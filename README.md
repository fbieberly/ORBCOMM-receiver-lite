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
    1. ```gcc -o ./decimator ./src/decimator.c -lm -lliquid -O2```
    1. ```gcc -o ./orbcomm_decoder ./src/orbcomm_decoder.c ./src/SGP4.c ./src/TLE.c ./src/eci2aer.c -lm -lliquid -O2```

## Real-time decoding
1. Update the orbcomm.txt TLE file from Celestrak.
    1. Run: ```wget -N https://www.celestrak.com/NORAD/elements/orbcomm.txt``` in the ORBCOMM-receiver-lite folder.
    1. I included an old TLE file. It will probably just work, but update when you can.
2. Update latitude, longitude, altitude and minimum horizon of your receiver in _config.h_
3. Run ```rtl_sdr -s 1200000 -g 0 -f 137524000 - | ./decimator | ./orbcomm_decoder```
    1. I recommend you use [gPredict] to know where the ORBCOMM satellites are.
    1. Note that not all the ORBCOMM satellites still transmit.
    1. All the decoded packets are printed to STDOUT. Pipe it to a file or consume it with another program if you want to save the data.  
  
[gPredict]: https://github.com/csete/gpredict

## Output  

Here is the output of the program:  
Note: the ### on the line means that packet does not pass the checksum OR has an unrecognized header.  
```
me@laptop:~/ORBCOMM-receiver-lite/src $ rtl_sdr -s 1200000 -g 0 -f 137524000 - | ./decimator | ./orbcomm_decoder 
Len of orbcomm_sats: 11
Found 1 device(s):
New Satellite FM116, Elevation:  6.9, Doppler: 1087.9 Hz
  0:  Realtek, RTL2838UHIDIR, SN: 00000001

Using device 0: Generic RTL2832U OEM
Found Rafael Micro R820T tuner
[R82XX] PLL not locked!
Sampling at 1200000 S/s.
Tuned to 137524000 Hz.
Tuner gain set to automatic.
Reading samples in async mode...
Allocating 15 zero-copy buffers
FM116, Chan: 0, 1eae5492d6697134638bd2aa
FM116, Chan: 0, 0a800dbc1c4042eb0222b947
FM116, Chan: 1, 0a80c4621d47172d0212d8bc
FM116, Chan: 0, 0a842778164c298b0222a6f3
FM116, Chan: 1, 0b01c610769899845302336b
FM116, Chan: 0, 0a888fc3073f32950122e309
FM116, Chan: 1, 0b013327a9bbbb884701dbd0
FM116, Chan: 0, 0a8ca7bc1f6c894f0112eaa7
FM116, Chan: 1, 0b01c8264355558cb1028258
FM116, Chan: 0, 0a907542233e82c50112975d
FM116, Chan: 1, 0b01152276888890b70039b7
FM116, Chan: 0, 0a94e9fc183d42dd01229e48
FM116, Chan: 1, 0b011629768888940f01dead
FM116, Chan: 0, 0a983f8b196b79ab0012b525
FM116, Chan: 1, 0b0168294355559887001c3b
FM116, Chan: 0, 0b0135276377779c1900dbb7
FM116, Chan: 1, 0b0169283133339c9700cfca
FM116, Chan: 0, 0aa0fdb1153ca23b002290c8
FM116, Chan: 1, 0b01771d999999a041015f54
FM116, Chan: 0, 0aa48c392146a7d502220a7c
FM116, Chan: 1, 0b013727768888a42d022914
FM116, Chan: 0, 0aa889051a3b325502127957
FM116, Chan: 1, 0aa8e1e5193af2e3021276d6
FM116, Chan: 0, 0aac0a70242b1ecb0212166e
FM116, Chan: 1, 0aac78b51c0757cf0132b6eb
FM116, Chan: 0, 0ab028651a24574f01127e44
FM116, Chan: 1, 0ab0b26325b94ee30012c947
FM116, Chan: 0, 0b015e29657777b4c501e6ba
FM116, Chan: 1, 1e567dc3683187453d728dab
FM116, Chan: 0, 0e01714d1002f01400007f9e
FM116, Chan: 1, 1e463e5d5fb36a15e6800109
FM116, Chan: 0, 0abcb5de7aeea6ab0022d5f7
FM116, Chan: 1, 1a1080d1820000000000ac57
FM116, Chan: 0, 1a104061012c6f080000a0f1
FM116, Chan: 1, 1a1011c5840000000000e597
FM116, Chan: 0, 0ac453831d1e18eb0222fa00
FM116, Chan: 1, 1a2062d31a762d30ad827203
FM116, Chan: 0, 1a10510100b137d80583b389
FM116, Chan: 1, 1a2142419500000000002b82
FM116, Chan: 0, 1a1092c33c33322898c24c12
FM116, Chan: 1, 1a10331d030b2d20ba061dce ### 
FM116, Chan: 0, 1a103328003a130000007cb2
FM116, Chan: 1, 1ec8eed0cd49cff8ebadb136
FM116, Chan: 0, 0e01714d1002f01400007f9e
FM116, Chan: 1, 1ea1f27df45a07ef620d1b04
```


## References

I used these two resources as my primary references.

http://mdkenny.customer.netspace.net.au/Orbcomm.pdf  
http://www.decodesystems.com/orbcomm.html  
