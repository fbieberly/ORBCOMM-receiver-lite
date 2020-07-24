// gcc -g -o ./decimator ./decimator.c -lm -lliquid

#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <liquid/liquid.h>

// Constants
#define SAMPLE_RATE 1.2e6
#define DECIM 2

int main(void)
{
    int chunk_size = 4096*32;
    uint8_t *data_buffer = (uint8_t *) malloc(chunk_size);

    // decimating LPF parameters
    float ft = 0.5;
    float As = 40.0f;
    float mu = 0.0f;
    unsigned int lpf_h_len = estimate_req_filter_len(ft, As);
    float lpf_h[lpf_h_len];

    // create taps for fir low pass filter
    liquid_firdes_kaiser(lpf_h_len, fc, As, mu, lpf_h);

    // create liquid decimating low pass filter object
    firdecim_crcf decim_lpf_q = firdecim_crcf_create(DECIM, lpf_h, lpf_h_len);
    firdecim_crcf_set_scale(decim_lpf_q, 1.0f/DECIM);

    // buffers for samples
    float complex samples[chunk_size/2];
    float complex decim_samples[chunk_size/4];

    while(1)
    {
        // Read samples in from STDIN
        int nsamps = read(STDIN_FILENO, data_buffer, chunk_size);
        if (nsamps < 10){ break; }

        // Convert uint8_t samples from RTLSDR to complex floats
        for (int j = 0; j < chunk_size/2; ++j)
        {
            float real = (float) data_buffer[2*j]/127.5 - 1.0;
            float imag = (float) data_buffer[2*j+1]/127.5 - 1.0;
            samples[j] = real + imag*_Complex_I;
        }

        // Low pass filter and decimate by 2
        firdecim_crcf_execute_block(decim_lpf_q, samples, chunk_size/2/DECIM, decim_samples);

        // write out samples to STDOUT
        fwrite(decim_samples, sizeof(decim_samples), 1, stdout);
    }

    free(data_buffer);
    return 0;
}
