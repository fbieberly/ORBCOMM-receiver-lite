# ORBCOMM receiver lite  

A software receiver for ORBCOMM satellite transmissions.  

## Description

This is a software receiver for decoding packets from ORBCOMM satellites. This program is 
a lightweight C implementation of my [ORBCOMM receiver]. I wrote it with the idea of 
being able to run it on a Raspberry Pi (3B or later). Additionally, this program will decode
all of the packets from both channels of all ORBCOMM satellites that are overhead.    

If you want a more full-featured ORBCOMM receiver please check out:  
https://www.coaa.co.uk/orbcommplotter.htm  
http://f6cte.free.fr/index_anglais.htm  

[ORBCOMM receiver]: https://github.com/fbieberly/ORBCOMM-receiver  

## Dependencies

The receiver takes samples from the rtl_sdr program that comes with [librtlsdr]. Additionally,
I use [liquid-dsp] for the signal processing.

[librtlsdr]: https://github.com/steve-m/librtlsdr
[liquid-dsp]: https://github.com/jgaeddert/liquid-dsp/

## Building

1. Install librtlsdr
2. Install liquid-dsp
3. Build this project
    1. ```gcc -o ./orbcomm_decoder ./orbcomm_decoder.c ./SGP4.c ./TLE.c ./eci2aer.c -lm -lliquid -O2```

## Real-time decoding
1. Update the orbcomm.txt TLE file from Celestrak.
    1. Run: ```wget -N https://www.celestrak.com/NORAD/elements/orbcomm.txt``` in the orbcomm_liquid folder.
2. Update latitude, longitude, altitude and minimum horizon of your receiver in _config.h_
3. Run _rtl_sdr -s 1228800 -g 0 -f 137500000 - | ./orbcomm_decoder_
    1. I recommend you use [gPredict] to know where the ORBCOMM satellites are.
    1. Note that not all the ORBCOMM satellites still transmit.
    1. All the decoded packets are printed to STDOUT. Pipe it to a file or consume it with another program if you want to save the data.  
  

## References

I used these two resources as my primary references.

http://mdkenny.customer.netspace.net.au/Orbcomm.pdf  
http://www.decodesystems.com/orbcomm.html  
