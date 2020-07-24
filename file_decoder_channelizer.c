// gcc -g -o ./decoder ./file_decoder_nodecim.c ./SGP4.c ./TLE.c ./eci2aer.c -lm -lliquid

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
#include "orbcomm_decoder.h"
#include "eci2aer.h"

// Constants
#define C 299792458.0
#define PI 3.14159
#define CENTER_FREQ 137500000.0
// #define SAMPLE_RATE 1.2288e6
#define SAMPLE_RATE 0.6144e6
#define DECIM 64

// Observer's location
// 43.802958, -99.210726, Pukwana
#define OBS_LAT 43.802958
#define OBS_LON -99.210726
#define OBS_ALT 0
#define MIN_ELEVATION 0

// # Derived from wikipedia: https://en.wikipedia.org/wiki/Fletcher%27s_checksum
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
    clock_t begin = clock();
    FILE *tle_file = NULL;
    FILE *data_file = NULL;
    FILE *out_file = NULL;
    char out_filename[] = "./out_samps.dat";

    double eci[3];
    double v[3];
    double aer[3];
    double lla[] = {OBS_LAT, OBS_LON, OBS_ALT};

    time_t rawtime;
    struct tm * ptm;
    double jd;
    double jdf;
    int hour_daylight_savings;

    // Set recording time and file name
    // rawtime = 1592155265.88; char samp_file[] = "./file1.dat";
    // rawtime = 1592155344.47; char samp_file[] = "./file2.dat";
    // rawtime = 1592155398.55; char samp_file[] = "./file3.dat";
    // rawtime = 1592155483.59; char samp_file[] = "./file4.dat";
    // rawtime = 1592155762.97; char samp_file[] = "./file5.dat";
    rawtime = 1593645769.09; char samp_file[] = "./file6.dat";

    ptm = gmtime(&rawtime);

    hour_daylight_savings = ptm->tm_hour;
    if (ptm->tm_isdst == 1){
        hour_daylight_savings--;
    }
    jday(ptm->tm_year+1900, ptm->tm_mon+1,
         ptm->tm_mday, hour_daylight_savings,
         ptm->tm_min, ptm->tm_sec,
         &jd, &jdf);

    int num_sats = sizeof(orbcomm_sats)/sizeof(orbcomm_sats[0]);
    printf("Len of orbcomm_sats: %ld\n", sizeof(orbcomm_sats)/sizeof(orbcomm_sats[0]));

    // Open TLE file and fill in data into struct
    tle_file = fopen("./orbcomm.txt","r");
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
            }
        }
    }
    fclose(tle_file);


    SAT_INFO sats_overhead[4];
    int num_sats_overhead = 0;

    int chunk_size = 4096*16;
    uint8_t *data_buffer = (uint8_t *) malloc(chunk_size*4);

    // Output samples to a file
    out_file = fopen(out_filename, "wb");

    // decimating LPF parameters
    float fc = 9600.0/SAMPLE_RATE;
    float ft = fc;
    float As = 60.0f;
    float mu = 0.0f;
    unsigned int lpf_h_len = estimate_req_filter_len(ft, As);
    float lpf_h[lpf_h_len];                   // filter coefficients
    liquid_firdes_kaiser(lpf_h_len, fc, As, mu, lpf_h);

    // firdecim_crcf decim_lpf_q = firdecim_crcf_create(DECIM, lpf_h, lpf_h_len);
    // firdecim_crcf_set_scale(decim_lpf_q, 1.0f/DECIM);

    // RRC pulse shape filter
    int samples_per_symbol = 2;
    int num_of_symbols_half_filter = 8;
    float rrc_alpha = 0.4;
    float bandwidth = 0.3;
    int rrc_num_taps = samples_per_symbol * num_of_symbols_half_filter * 2 + 1;

    // symtrack_cccf symtrack_q = symtrack_cccf_create(LIQUID_FIRFILT_RRC,
    //                                                samples_per_symbol,
    //                                                num_of_symbols_half_filter,
    //                                                rrc_alpha,
    //                                                LIQUID_MODEM_BPSK);
    // symtrack_cccf_set_bandwidth(symtrack_q, bandwidth);


    // symsync_crcf symsync_q = symsync_crcf_create_rnyquist(LIQUID_FIRFILT_RRC,
    //                                                       samples_per_symbol,
    //                                                       num_of_symbols_half_filter,
    //                                                       rrc_alpha,
    //                                                       31);


    int loop_counter = 0;
    // int chunk_size = 4096*8;
    int num_symbols = chunk_size/2;
    int num_bits = num_symbols/samples_per_symbol;
    int num_syms = 0;

    float complex last_symbol = 1 + 1*_Complex_I;

    uint8_t bits[2000];
    int bits_start = 0;
    int bits_stop = 0;

    uint8_t num_stages = 7;
    int r = 1<<num_stages;

    float complex samples[chunk_size/2];


    while(1)
    {
        for (int j = 0; j < num_sats; ++j)
        {
            TLE tle;
            double eci[3];
            double v[3];
            double aer[3];
            parseLines(&tle, orbcomm_sats[j].line1, orbcomm_sats[j].line2);
            getRVForDate(&tle, rawtime*1000.0, eci, v);
            eci2aer(jd+jdf, eci, lla, aer);

            if (aer[1] >= MIN_ELEVATION){
                // calculate doppler by calculating the range to the
                // satellite 1 second in the future and differencing
                // the slant ranges.
                double old_range = aer[2];
                getRVForDate(&tle, rawtime*1000.0 + 1000.0, eci, v);
                eci2aer(jd+jdf+1.0/(24*60.0*60.0), eci, lla, aer);
                double delta_range =  aer[2] - old_range;
                double doppler = C/(C+delta_range) * 137500000.0 - 137500000.0;
                // printf("AER: %f %f %f\n", aer[0], aer[1], aer[2]);
                // printf("Relative range rate: %6.1f m/s\n", delta_range);

                orbcomm_sats[j].doppler = doppler;
                uint8_t found_sat = 0;
                for (int k = 0; k < num_sats_overhead; ++k)
                {
                    char * ind = strstr(sats_overhead[k].sat_name, orbcomm_sats[j].sat_name);
                    if (ind != NULL)
                    {
                        sats_overhead[k].doppler = doppler;
                        found_sat = 1;
                        // printf("Update Satellite %s, Elevation: %4.1f, Doppler: %6.1f Hz\n", orbcomm_sats[j].sat_name, orbcomm_sats[j].aer[1], doppler);
                        break;
                    }
                }
                if (found_sat == 0){
                    sats_overhead[num_sats_overhead] = orbcomm_sats[j];
                    printf("New Satellite %s, Elevation: %4.1f, Doppler: %6.1f Hz\n", orbcomm_sats[j].sat_name, aer[1], doppler);
                    // sats_overhead[num_sats_overhead].chan1_decoder.nco_q = nco_crcf_create(LIQUID_VCO);
                    // nco_crcf_set_phase(sats_overhead[num_sats_overhead].chan1_decoder.nco_q, 0.0f);
                    // nco_crcf_set_frequency(sats_overhead[num_sats_overhead].chan1_decoder.nco_q, 0.0f);

                    sats_overhead[num_sats_overhead].decoder[0].nco_q = nco_crcf_create(LIQUID_VCO);
                    nco_crcf_set_phase(sats_overhead[num_sats_overhead].decoder[0].nco_q, 0.0f);
                    nco_crcf_set_frequency(sats_overhead[num_sats_overhead].decoder[0].nco_q, 0.0f);

                    sats_overhead[num_sats_overhead].decoder[1].nco_q = nco_crcf_create(LIQUID_VCO);
                    nco_crcf_set_phase(sats_overhead[num_sats_overhead].decoder[1].nco_q, 0.0f);
                    nco_crcf_set_frequency(sats_overhead[num_sats_overhead].decoder[1].nco_q, 0.0f);

                    sats_overhead[num_sats_overhead].decoder[0].decim_lpf_q = firdecim_crcf_create(DECIM, lpf_h, lpf_h_len);
                    firdecim_crcf_set_scale(sats_overhead[num_sats_overhead].decoder[0].decim_lpf_q, 1.0f/DECIM);

                    sats_overhead[num_sats_overhead].decoder[1].decim_lpf_q = firdecim_crcf_create(DECIM, lpf_h, lpf_h_len);
                    firdecim_crcf_set_scale(sats_overhead[num_sats_overhead].decoder[1].decim_lpf_q, 1.0f/DECIM);

                    sats_overhead[num_sats_overhead].decoder[0].symsync_q = symsync_crcf_create_rnyquist(LIQUID_FIRFILT_RRC,
                                                          samples_per_symbol,
                                                          num_of_symbols_half_filter,
                                                          rrc_alpha,
                                                          31);
                    sats_overhead[num_sats_overhead].decoder[1].symsync_q = symsync_crcf_create_rnyquist(LIQUID_FIRFILT_RRC,
                                                          samples_per_symbol,
                                                          num_of_symbols_half_filter,
                                                          rrc_alpha,
                                                          31);

                    // printf("chunk/128: %f\n", (float)chunk_size/2/r);
                    float shift_freq_hz = CENTER_FREQ - sats_overhead[num_sats_overhead].channel[0] - sats_overhead[num_sats_overhead].doppler;
                    float shift_fc = shift_freq_hz/SAMPLE_RATE;
                    sats_overhead[num_sats_overhead].decoder[0].dds_cccf_q = dds_cccf_create(num_stages, -shift_fc, 0.25, 60.0);

                    shift_freq_hz = CENTER_FREQ - sats_overhead[num_sats_overhead].channel[1] - sats_overhead[num_sats_overhead].doppler;
                    shift_fc = shift_freq_hz/SAMPLE_RATE;
                    sats_overhead[num_sats_overhead].decoder[1].dds_cccf_q = dds_cccf_create(num_stages, -shift_fc, 0.25, 60.0);

                    num_sats_overhead++;

                }
            }
        }

        // Read samples in from STDIN
        int nsamps = read(STDIN_FILENO, data_buffer, chunk_size*4);
        if (nsamps < 100){ break; }

        float * buff_ptr = (float *) data_buffer;
        for (int j = 0; j < chunk_size/2; ++j)
        {
            // float real = (float) data_buffer[2*j]/127.5 - 1.0;
            // float imag = (float) data_buffer[2*j+1]/127.5 - 1.0;
            float real = (float) buff_ptr[j];
            float imag = (float) buff_ptr[j+1];
            buff_ptr += 1;
            samples[j] = real + imag*_Complex_I;
        }


        fwrite(samples, sizeof(samples), 1, out_file);

        // continue;


        /* Do stuff with samples */
        // normalize samples by the median(abs(samples))
        // float complex agc_samples[chunk_size/2];
        // for (int j = 0; j < chunk_size/2; ++j)
        // {
        //     agc_crcf_execute(agc_q, samples[j], &agc_samples[j]);
        // }

        float complex bb_samples[chunk_size/2];
        float complex bb_decim_samples[chunk_size/2/DECIM];
        float complex bb_decim_samples_2[chunk_size/2/r];
        float complex rrc_samples[chunk_size/2/DECIM/2];

        float angles[chunk_size/2/DECIM/2];

        unsigned int num_written;


        num_sats_overhead = 1;
        uint8_t channel_num = 0;


        for (int satnum = 0; satnum < num_sats_overhead; ++satnum)
        {
            // frequency shift to baseband, use updated doppler shift
            float shift_freq_hz = CENTER_FREQ - sats_overhead[satnum].channel[channel_num] - sats_overhead[satnum].doppler;
            // nco_crcf_set_frequency(nco_q, shift_freq_hz/SAMPLE_RATE*2*PI);
            // nco_crcf_mix_block_up(nco_q, agc_samples, bb_samples, chunk_size/2);
            nco_crcf_set_frequency(sats_overhead[satnum].decoder[channel_num].nco_q, shift_freq_hz/SAMPLE_RATE*2*PI);
            nco_crcf_mix_block_up(sats_overhead[satnum].decoder[channel_num].nco_q, samples, bb_samples, chunk_size/2);

            // // decimate by 128/low pass filter
            firdecim_crcf_execute_block(sats_overhead[satnum].decoder[channel_num].decim_lpf_q, bb_samples, chunk_size/2/DECIM, bb_decim_samples);

            // printf("Frequency: %f\n", shift_fc);
            // printf("r: %d\n", r);

            // dds_cccf_print(dds_cccf_q);

            // run decimation (down-conversion) stage
            // for (int j = 0; j < chunk_size/2/r; j++) {
            //     dds_cccf_decim_execute(sats_overhead[satnum].decoder[channel_num].dds_cccf_q, &agc_samples[r*j], &bb_decim_samples_2[j]);
            // }

            // RRC matched filter, timing recovery
            // All in one!!
            symsync_crcf_execute(sats_overhead[satnum].decoder[channel_num].symsync_q, bb_decim_samples, chunk_size/2/DECIM,
                                 rrc_samples, &num_written);
            // printf("Num Written: %d\n", num_written);
            // printf("Num Written: %d\n", chunk_size/2/DECIM);
            // Write some samples to disk so we can plot
            // fwrite(bb_decim_samples_2, sizeof(bb_decim_samples_2), 1, out_file);

            loop_counter++;
            num_syms += (chunk_size/2/DECIM/2);

            // demod symbols
            float tmp_angle;
            for (int j = 0; j < chunk_size/2/DECIM/2; ++j)
            {
                tmp_angle = (carg(rrc_samples[j]) - carg(last_symbol)) * 180.0 / PI;
                if (tmp_angle >  180.0){ tmp_angle -= 360.0; }
                if (tmp_angle < -180.0){ tmp_angle += 360.0; }
                angles[j] = tmp_angle;
                last_symbol = rrc_samples[j];
            }

            // produce bits
            uint8_t bit = 0;
            for (int j = 0; j < chunk_size/2/DECIM/2; ++j)
            {
                if (angles[j] > 0){ bit = 1; }
                else { bit = 0; }
                bits[bits_stop + j] = bit;
            }
            bits_stop += chunk_size/2/DECIM/2;
            // printf("Bit stop: %d\n", bits_stop);



            uint8_t packet[12];
            uint8_t header_found = 0;
            uint8_t last_packet_bad = 0;
            while (bits_stop - bits_start > 96){

                for (int k = 0; k < 12; ++k)
                {
                    uint8_t packbits = 0;
                    packbits += bits[k*8+0+bits_start];
                    packbits += bits[k*8+1+bits_start] << 1;
                    packbits += bits[k*8+2+bits_start] << 2;
                    packbits += bits[k*8+3+bits_start] << 3;
                    packbits += bits[k*8+4+bits_start] << 4;
                    packbits += bits[k*8+5+bits_start] << 5;
                    packbits += bits[k*8+6+bits_start] << 6;
                    packbits += bits[k*8+7+bits_start] << 7;
                    packet[k] = packbits;
                }
                header_found = 0;
                for (int k = 0; k < sizeof(packet_headers)/sizeof(packet_headers[0]); ++k){
                    if ( packet[0] == packet_headers[k]){
                        header_found = 1;
                        break;
                    }
                }

                uint16_t chcksum = fletcher16(packet, 12);
                if (chcksum == 0 && header_found == 1){
                    last_packet_bad = 0;

                    for (int k = 0; k < 12; ++k){
                        printf("%02x", packet[k]);
                    }
                    printf("\n");

                    bits_start += 96;

                } else { bits_start++; }
            }
            if (bits_start > 0){
                for (int j = 0; j < (bits_stop - bits_start); ++j){
                    bits[j] = bits[bits_start + j];
                }
                bits_stop -= bits_start;
                bits_start = 0;
            }


            // fwrite(bits, sizeof(bits), 1, out_file);

            break;
        }
    }
    printf("\n");
    printf("loop_counter: %d\n", loop_counter);
    // printf("fill count: %d\n", fill_count);

    // agc_crcf_destroy(agc_q);
    // nco_crcf_destroy(nco_q);
    // firdecim_crcf_destroy(decim_lpf_q);
    // symtrack_cccf_destroy(symtrack_q);
    // symsync_crcf_destroy(symsync_q);


    fclose(out_file);
    free(data_buffer);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Time taken: %f\n", time_spent);
    return 0;
}
