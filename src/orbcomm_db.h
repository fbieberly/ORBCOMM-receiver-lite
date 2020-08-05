

// Two decoders per satellite. One for each channel.
typedef struct ORBCOMM_DECODER{
    nco_crcf nco_q;
    symtrack_cccf symtrack_q;
    symsync_crcf symsync_q;
    firdecim_cccf firdecim_q;
    uint8_t bits[2000];
    int bits_start;
    int bits_stop;
    float complex last_symbol;
    uint8_t bad_packet_count;
    uint8_t offset_search;
} ORBCOMM_DECODER;


// Contains satellite info and decoders.
typedef struct SAT_INFO {
    char sat_name[255];
    char line1[255];
    char line2[255];
    float channel[2];
    float offset[2];
    float doppler;
    ORBCOMM_DECODER decoder[2];
    uint8_t overhead;
} SAT_INFO;

SAT_INFO orbcomm_sats[] = {
    { .sat_name = "FM103",
     .channel[0] = 137250000.0,
     .channel[1] = 137312500.0
    },

    { .sat_name = "FM107",
     .channel[0] = 137250000.0,
     .channel[1] = 137312500.0
    },

    { .sat_name = "FM108",
     .channel[0] = 137460000.0,
     .channel[1] = 137712500.0
    },

    { .sat_name = "FM109",
     .channel[0] = 137250000.0,
     .channel[1] = 137312500.0
    },

    { .sat_name = "FM110",
     .channel[0] = 137287500.0,
     .channel[1] = 137737500.0
    },

    { .sat_name = "FM112",
     .channel[0] = 137662500.0,
     .channel[1] = 137800000.0
    },

    { .sat_name = "FM113",
     .channel[0] = 137662500.0,
     .channel[1] = 137800000.0
    },

    { .sat_name = "FM114",
     .channel[0] = 137287500.0,
     .channel[1] = 137737500.0
    },

    { .sat_name = "FM116",
     .channel[0] = 137662500.0,
     .channel[1] = 137800000.0
    },

    { .sat_name = "FM117",
     .channel[0] = 137460000.0,
     .channel[1] = 137712500.0
    },

    { .sat_name = "FM118",
     .channel[0] = 137287500.0,
     .channel[1] = 137737500.0
    },
};

int packet_headers[] = {
    0x0a, // unkn 1
    0x0b, // unkn 2
    0x0c, // unkn 3
    0x0d, // unkn 4
    0x0e, // unkn 5
    0x0f, // unkn 6
    0x13, // unkn 7
    0x18, // unkn 8
    0x1a, // msg
    0x1b, // uplink
    0x1c, // downlink
    0x1d, // network
    0x1e, // fill
    0x1f, // ephemeris
    0x22, // orbital
    0x65, // sync
};