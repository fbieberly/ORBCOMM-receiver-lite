// gcc -I -o ./orbcomm_decoder ./src/orbcomm_decoder.c ./src/SGP4.c ./src/TLE.c ./src/eci2aer.c -lm -lliquid -O2

#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// headers from someone else
#include <liquid/liquid.h>
#include "TLE.h"

// Files I made
#include "orbcomm_db.h"
#include "eci2aer.h"
#include "../config.h"

// Constants
#define C 299792458.0
#define CENTER_FREQ 137524000.0
#define SAMPLE_RATE 0.6e6
#define DECIM 25

// Derived from wikipedia: https://en.wikipedia.org/wiki/Fletcher%27s_checksum
// Note: This implementation is different than the classical fletcher checksum
//       as it uses a 256 modulo instead of 255.
uint16_t fletcher16( uint8_t *data, int count ){
   uint16_t sum1 = 0;
   uint16_t sum2 = 0;
   int index;
   for ( index = 0; index < count; ++index ){
      sum1 = (sum1 + data[index]) % 256;
      sum2 = (sum2 + sum1) % 256;
   }
   return (sum2 << 8) | sum1;
}

int main(void)
{
    FILE *tle_file = NULL;
    double lla[] = {OBS_LAT, OBS_LON, OBS_ALT};
    time_t rawtime;
    time_t old_rawtime = 0.0;;
    double time_increase = 0.0;

    struct tm * ptm;
    double jd, jdf;
    int hour_daylight_savings;

    int num_sats = sizeof(orbcomm_sats)/sizeof(orbcomm_sats[0]);
    printf("Len of orbcomm_sats: %d\n", num_sats);

    // Open TLE file and fill in data into struct
    tle_file = fopen("../orbcomm.txt","r");
    for (int i = 0; i < num_sats; ++i)
    {
        char sat_name[255];
        strcpy(sat_name, orbcomm_sats[i].sat_name);
        fseek( tle_file, 0, SEEK_SET );

        char line0[255];
        while(fgets(line0, 255, tle_file) != NULL)
        {
            char * ind = strstr(line0, sat_name);
            if (ind != NULL)
            {
                ind = fgets(orbcomm_sats[i].line1, 255, tle_file);
                ind = fgets(orbcomm_sats[i].line2, 255, tle_file);
                break;
            }
        }
    }
    fclose(tle_file);

    int chunk_size = 250 * 48;
    uint8_t *data_buffer = (uint8_t *) malloc(chunk_size * sizeof(float));

    // decimating LPF parameters
    float fc = 19200.0/SAMPLE_RATE;
    float As = 40.0f;
    float mu = 0.0f;
    unsigned int lpf_h_len = estimate_req_filter_len(fc, As);

    float lpf_h[lpf_h_len];                   // filter coefficients
    liquid_firdes_kaiser(lpf_h_len, fc, As, mu, lpf_h);
    float complex bb_bpf_h[lpf_h_len];
    float complex mixed_bpf_h[lpf_h_len];
    for (int i = 0; i < lpf_h_len; ++i)
    {
        bb_bpf_h[i] = lpf_h[i] + lpf_h[i]*_Complex_I;
    }
    nco_crcf filt_nco_q = nco_crcf_create(LIQUID_VCO);

    // RRC pulse shape filter parameters
    int samples_per_symbol = 5;
    int num_of_symbols_half_filter = 8;
    float rrc_alpha = 0.4;
    float bandwidth = 0.3;
    int rrc_num_taps = samples_per_symbol * num_of_symbols_half_filter * 2 + 1;

    // Initialize liquid objects for each decoder
    for (int i = 0; i < num_sats; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            // NCO
            orbcomm_sats[i].decoder[j].nco_q = nco_crcf_create(LIQUID_VCO);
            nco_crcf_reset(orbcomm_sats[i].decoder[j].nco_q);

            // Decimating BPF
            float freq_shift = (orbcomm_sats[i].channel[j] - CENTER_FREQ)/SAMPLE_RATE;
            nco_crcf_set_phase(filt_nco_q, 0.0f);
            nco_crcf_set_frequency(filt_nco_q, freq_shift*2*M_PI);
            nco_crcf_mix_block_up(filt_nco_q, bb_bpf_h, mixed_bpf_h, lpf_h_len);
            orbcomm_sats[i].decoder[j].firdecim_q = firdecim_cccf_create(DECIM, mixed_bpf_h, lpf_h_len);

            // Symbol tracking
            orbcomm_sats[i].decoder[j].symtrack_q = symtrack_cccf_create(LIQUID_FIRFILT_RRC,
                                                                         samples_per_symbol,
                                                                         num_of_symbols_half_filter,
                                                                         rrc_alpha,
                                                                         LIQUID_MODEM_QPSK);
        }
    }
    nco_crcf_destroy(filt_nco_q);

    float complex samples[2*chunk_size];
    float complex filtered_samples[2*chunk_size];
    float complex bb_decim_samples[chunk_size/2/DECIM];
    float complex rrc_samples[chunk_size/2/DECIM/samples_per_symbol];
    float angles[chunk_size/2/DECIM/samples_per_symbol];
    unsigned int num_written;

    TLE tle;
    double eci[3], v[3], aer[3];

    float decim_top = CENTER_FREQ + SAMPLE_RATE/DECIM/2;
    float decim_bot = CENTER_FREQ - SAMPLE_RATE/DECIM/2;

    int sample_start = 0;
    int sample_stop = 0;
    int total_samples = 0;
    uint8_t found_sat = 0;
    clock_t begin = clock();

    while (1)
    {
        time ( &rawtime );

        if (rawtime != old_rawtime){
            old_rawtime = rawtime;
            time_increase = 0.0;
        }
        ptm = gmtime(&rawtime);

        hour_daylight_savings = ptm->tm_hour;
        if (ptm->tm_isdst == 1){
            hour_daylight_savings--;
        }
        jday(ptm->tm_year+1900, ptm->tm_mon+1,
             ptm->tm_mday, hour_daylight_savings,
             ptm->tm_min, ptm->tm_sec+time_increase,
             &jd, &jdf);

        found_sat = 0;
        for (int j = 0; j < num_sats; ++j)
        {
            parseLines(&tle, orbcomm_sats[j].line1, orbcomm_sats[j].line2);
            getRVForDate(&tle, (rawtime + time_increase)*1000.0, eci, v);
            eci2aer(jd+jdf, eci, lla, aer);

            if (aer[1] >= MIN_ELEVATION){
                found_sat = 1;

                // calculate doppler by calculating the range to the satellite 1 second in
                // the future and differencing the slant ranges.
                double old_range = aer[2];
                getRVForDate(&tle, (rawtime + time_increase)*1000.0 + 1000.0, eci, v);
                eci2aer(jd+jdf+1.0/(24.0*60.0*60.0), eci, lla, aer);
                double delta_range =  aer[2] - old_range;
                double doppler = C/(C+delta_range) * CENTER_FREQ - CENTER_FREQ;
                orbcomm_sats[j].doppler = doppler;
                // printf("AER: %f %f %f\n", aer[0], aer[1], aer[2]);
                // printf("Relative range rate: %6.1f m/s\n", delta_range);

                if (orbcomm_sats[j].overhead == 0){
                    printf("New Satellite %s, Elevation: %4.1f, Doppler: %6.1f Hz\n", orbcomm_sats[j].sat_name, aer[1], doppler);
                    orbcomm_sats[j].overhead = 1;

                    for (int k = 0; k < 2; ++k)
                    {
                        nco_crcf_reset(orbcomm_sats[j].decoder[k].nco_q);
                        symtrack_cccf_reset(orbcomm_sats[j].decoder[k].symtrack_q);

                        orbcomm_sats[j].decoder[k].bits_start = 0;
                        orbcomm_sats[j].decoder[k].bits_stop = 0;
                        orbcomm_sats[j].decoder[k].last_symbol = 1 + 1*_Complex_I;
                        orbcomm_sats[j].decoder[k].bad_packet_count = 0;
                        orbcomm_sats[j].decoder[k].offset_search = 1;
                    }
                }
            } else {
                if (orbcomm_sats[j].overhead == 1){
                    // printf("Satellite %s now below horizon.\n", orbcomm_sats[j].sat_name);
                    // for (int k = 0; k < 2; ++k){
                    //     nco_crcf_destroy(orbcomm_sats[j].decoder[k].nco_q);
                    //     firdecim_cccf_destroy(orbcomm_sats[j].decoder[k].firdecim_q);
                    //     symtrack_cccf_destroy(orbcomm_sats[j].decoder[k].symtrack_q);
                    // }
                }
                orbcomm_sats[j].overhead = 0;
            }
        }

        // Read samples in from STDIN
        int num_samps = read(STDIN_FILENO, data_buffer, chunk_size * sizeof(float));
        // if (num_samps <= 10) { break; }
        time_increase += (float) num_samps / SAMPLE_RATE;

        float * buff_ptr = (float *) data_buffer;
        for (int j = 0; j < num_samps/4; ++j)
        {
            float real = (float) buff_ptr[2*j];
            float imag = (float) buff_ptr[2*j+1];
            samples[sample_stop + j] = real + imag*_Complex_I;
        }
        sample_stop += num_samps/8;

        // When too few samples are read from stdin, continue
        if (sample_stop - sample_start < chunk_size/2) { continue; }

        // Else, we will process chunk_size/2 samples
        total_samples += chunk_size/2;

        for (int satnum = 0; satnum < num_sats; ++satnum)
        {
            if (orbcomm_sats[satnum].overhead)
            {
                for (int channel_num = 0; channel_num < 2; ++channel_num)
                {
                    // frequency shift to baseband, use updated doppler shift
                    float shift_freq_hz;
                    if (orbcomm_sats[satnum].channel[channel_num] < CENTER_FREQ){
                        float freq_below = decim_bot - orbcomm_sats[satnum].channel[channel_num] - orbcomm_sats[satnum].doppler;
                        float freq_offset = fmod(freq_below, SAMPLE_RATE/DECIM);
                        shift_freq_hz = CENTER_FREQ - (decim_top - freq_offset);
                    } else {
                        float freq_above = orbcomm_sats[satnum].channel[channel_num] + orbcomm_sats[satnum].doppler - decim_top;
                        float freq_offset = fmod(freq_above, SAMPLE_RATE/DECIM);
                        shift_freq_hz = CENTER_FREQ - (decim_bot + freq_offset);
                    }

                    firdecim_cccf_execute_block(orbcomm_sats[satnum].decoder[channel_num].firdecim_q, samples, chunk_size/2/DECIM, filtered_samples);
                    nco_crcf_set_frequency(orbcomm_sats[satnum].decoder[channel_num].nco_q,
                                           shift_freq_hz/(SAMPLE_RATE/DECIM)*2*M_PI);

                    nco_crcf_mix_block_up(orbcomm_sats[satnum].decoder[channel_num].nco_q,
                                          filtered_samples, bb_decim_samples, chunk_size/2/DECIM);

                    symtrack_cccf_execute_block(orbcomm_sats[satnum].decoder[channel_num].symtrack_q,
                                                bb_decim_samples,
                                                chunk_size/2/DECIM,
                                                rrc_samples,
                                                &num_written);

                    // demod symbols
                    float tmp_angle;
                    for (int j = 0; j < chunk_size/2/DECIM/samples_per_symbol; ++j){
                        tmp_angle = (carg(rrc_samples[j]) - carg(orbcomm_sats[satnum].decoder[channel_num].last_symbol)) * 180.0 / M_PI;
                        if (tmp_angle >  180.0){ tmp_angle -= 360.0; }
                        if (tmp_angle < -180.0){ tmp_angle += 360.0; }
                        angles[j] = tmp_angle;
                        orbcomm_sats[satnum].decoder[channel_num].last_symbol = rrc_samples[j];
                    }

                    int bits_start = orbcomm_sats[satnum].decoder[channel_num].bits_start;
                    int bits_stop = orbcomm_sats[satnum].decoder[channel_num].bits_stop;

                    // produce bits
                    uint8_t bit = 0;
                    for (int j = 0; j < chunk_size/2/DECIM/samples_per_symbol; ++j)
                    {
                        if (angles[j] > 0){ bit = 1; }
                        else { bit = 0; }
                        orbcomm_sats[satnum].decoder[channel_num].bits[bits_stop + j] = bit;
                    }
                    bits_stop += chunk_size/2/DECIM/samples_per_symbol;

                    uint8_t packet_len = 12;
                    uint8_t header_found;
                    uint8_t packet[24];
                    while (bits_stop - bits_start > 24*8){
                        packet_len = 12;
                        // create a 12 byte packet
                        for (int k = 0; k < packet_len; ++k){
                            uint8_t packbits = 0;
                            packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+0+bits_start];
                            packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+1+bits_start] << 1;
                            packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+2+bits_start] << 2;
                            packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+3+bits_start] << 3;
                            packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+4+bits_start] << 4;
                            packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+5+bits_start] << 5;
                            packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+6+bits_start] << 6;
                            packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+7+bits_start] << 7;
                            packet[k] = packbits;
                        }
                        // check for a known header
                        header_found = 0;
                        for (int k = 0; k < sizeof(packet_headers)/sizeof(packet_headers[0]); ++k){
                            if ( packet[0] == packet_headers[k]){
                                header_found = 1;
                                break;
                            }
                        }
                        // For certain packets (e.g. ephemeris) the packet is 24 bytes long
                        if (packet[0] == 0x1f || packet[0] == 0x21){
                            packet_len = 24;
                            for (int k = 0; k < packet_len; ++k){
                                uint8_t packbits = 0;
                                packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+0+bits_start];
                                packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+1+bits_start] << 1;
                                packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+2+bits_start] << 2;
                                packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+3+bits_start] << 3;
                                packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+4+bits_start] << 4;
                                packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+5+bits_start] << 5;
                                packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+6+bits_start] << 6;
                                packbits += orbcomm_sats[satnum].decoder[channel_num].bits[k*8+7+bits_start] << 7;
                                packet[k] = packbits;
                            }
                            header_found = 1;
                        }

                        uint16_t checksum = fletcher16(packet, packet_len);
                        if (checksum == 0 && header_found == 1){
                            orbcomm_sats[satnum].decoder[channel_num].offset_search = 0;
                            orbcomm_sats[satnum].decoder[channel_num].bad_packet_count = 0;

                            printf("%s, ", orbcomm_sats[satnum].sat_name);
                            printf("Chan: %d, ", channel_num);
                            for (int k = 0; k < packet_len; ++k){
                                printf("%02x", packet[k]);
                            }
                            printf("\n");
                            bits_start += packet_len*8;
                        } else {
                            if (orbcomm_sats[satnum].decoder[channel_num].offset_search == 0){
                                orbcomm_sats[satnum].decoder[channel_num].bad_packet_count += 1;
                                bits_start += packet_len*8;
                                printf("%s, ", orbcomm_sats[satnum].sat_name);
                                printf("Chan: %d, ", channel_num);
                                for (int k = 0; k < packet_len; ++k){
                                    printf("%02x", packet[k]);
                                }
                                printf(" ### \n");
                            }
                        }

                        if (orbcomm_sats[satnum].decoder[channel_num].offset_search == 1) { bits_start++; }
                        if (orbcomm_sats[satnum].decoder[channel_num].offset_search == 0 && orbcomm_sats[satnum].decoder[channel_num].bad_packet_count >= 4){
                            orbcomm_sats[satnum].decoder[channel_num].offset_search = 1;
                        }
                    }
                    if (bits_start > 0){
                        for (int j = 0; j < (bits_stop - bits_start); ++j){
                            orbcomm_sats[satnum].decoder[channel_num].bits[j] = orbcomm_sats[satnum].decoder[channel_num].bits[bits_start + j];
                        }
                        bits_stop -= bits_start;
                        bits_start = 0;
                    }
                    orbcomm_sats[satnum].decoder[channel_num].bits_start = bits_start;
                    orbcomm_sats[satnum].decoder[channel_num].bits_stop = bits_stop;
                }
            }
        }
        sample_start += chunk_size/2;
        for (int j = 0; j < (sample_stop - sample_start); ++j){
            samples[j] = samples[sample_start + j];
        }
        sample_stop = sample_stop - sample_start;
        sample_start = 0;

        if (found_sat == 0){
            sleep(1);
        }
    }
    printf("\n");

    for (int j = 0; j < num_sats; ++j)
    {
        for (int k = 0; k < 2; ++k)
        {
            nco_crcf_destroy(orbcomm_sats[j].decoder[k].nco_q);
            firdecim_cccf_destroy(orbcomm_sats[j].decoder[k].firdecim_q);
            symtrack_cccf_destroy(orbcomm_sats[j].decoder[k].symtrack_q);
        }
    }

    free(data_buffer);
    // clock_t end = clock();
    // double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    // printf("Time per sample: %.20lf\n", time_spent/(total_samples*2));
    printf("Exiting.\n");
    return 0;
}
