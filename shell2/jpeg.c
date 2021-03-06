
#include <pogo.h>
#include "misc.h"
#include "iwram.h"
#include "viewer.h"

/* These macros are so that we can generate the AA&N multipliers at
 * compile-time, allowing configuration control of fixed point precision.
 */
#define JPEG_AAN_0 1.0
#define JPEG_AAN_1 1.387039845 
#define JPEG_AAN_2 1.306562965
#define JPEG_AAN_3 1.175875602
#define JPEG_AAN_4 1.0
#define JPEG_AAN_5 0.785694958
#define JPEG_AAN_6 0.541196100
#define JPEG_AAN_7 0.275899379

#define JPEG_AAN_LINE(B) \
    JPEG_FTOFIX (JPEG_AAN_0 * JPEG_AAN_##B), \
    JPEG_FTOFIX (JPEG_AAN_1 * JPEG_AAN_##B), \
    JPEG_FTOFIX (JPEG_AAN_2 * JPEG_AAN_##B), \
    JPEG_FTOFIX (JPEG_AAN_3 * JPEG_AAN_##B), \
    JPEG_FTOFIX (JPEG_AAN_4 * JPEG_AAN_##B), \
    JPEG_FTOFIX (JPEG_AAN_5 * JPEG_AAN_##B), \
    JPEG_FTOFIX (JPEG_AAN_6 * JPEG_AAN_##B), \
    JPEG_FTOFIX (JPEG_AAN_7 * JPEG_AAN_##B)

int jpegWidth, jpegHeight;

/* The AA&N scaling factors.  These should be multiplied against quantization
 * coefficients to determine their real value.
 */
const JPEG_FIXED_TYPE JPEG_AANScaleFactor [JPEG_DCTSIZE2] =
{
    JPEG_AAN_LINE (0),
    JPEG_AAN_LINE (1),
    JPEG_AAN_LINE (2),
    JPEG_AAN_LINE (3),
    JPEG_AAN_LINE (4),
    JPEG_AAN_LINE (5),
    JPEG_AAN_LINE (6),
    JPEG_AAN_LINE (7),
};


/* This function comes from jpeglib.  I feel all right about that since it comes from AA&N anyway. */
void JPEG_IDCT (JPEG_FIXED_TYPE *zz, signed char *chunk, int chunkStride)
{
    IDCT_Columns (zz);
    IDCT_Rows (zz, chunk, chunkStride);
}

/* Compute a signed value.  COUNT is the number of bits to read, and OUT is
 * where to store the result.
 */

#if 0 
#define JPEG_Value(COUNT, OUT) \
    do { \
        unsigned int value = JPEG_BITS_GET (COUNT); \
        \
        if (value < (unsigned int) (1 << ((unsigned int) (COUNT - 1)))) \
            value += (-1 << COUNT) + 1; \
        (OUT) = value; \
    } while (0)
#endif

/* Start decoding bits. */
unsigned int bits_left;
unsigned long int bits_data;

int bytes_left;
char *rewind_point;

#define decoder_value (0x02040000-sizeof(JPEG_Decoder))
JPEG_Decoder * decoder;
// = (JPEG_Decoder * const) (decoder_value); // sizeof(JPEG_Decoder)

/* The decompressed AC Huffman tables.  JPEG Baseline allows only two AC Huffman tables in a scan. */
#define acTableList_value (decoder_value-sizeof(JPEG_HuffmanTable)*2)
JPEG_HuffmanTable * acTableList;
// = (JPEG_HuffmanTable * const) (acTableList_value); // sizeof(JPEG_HuffmanTable)*2

/* The decompressed DC Huffman tables.  JPEG Baseline allows only two DC Huffman tables in a scan. */
#define dcTableList_value (acTableList_value-sizeof(JPEG_HuffmanTable)*2)
JPEG_HuffmanTable * dcTableList;
// = (JPEG_HuffmanTable * const) (dcTableList_value); // sizeof(JPEG_HuffmanTable)*2

#define jpegHeap (dcTableList_value)

void JPEG_Decoder_ReadImage_Init(const char *dataBase)
{
    JPEG_FrameHeader_Component *item, *itemEnd; /* The frame header's components for loops. */
    int horzShift = 0; /* The right shift to use after multiplying by nHorzFactor to get the actual sample. */
    int vertShift = 0; /* The right shift to use after multiplying by nVertFactor to get the actual sample. */
	int c;

	decoder->YHorzFactor = decoder->YVertFactor = 0;
    decoder->CbHorzFactor = decoder->CbVertFactor = decoder->CrHorzFactor = decoder->CrVertFactor = 1;
	decoder->horzMax = decoder->vertMax = 0;
    decoder->acTableUse[0] = decoder->acTableUse[1] = -1;
    decoder->dcTableUse[0] = decoder->dcTableUse[1] = -1;
    decoder->factorSum = 0;
	decoder->startdata = decoder->data = dataBase;
    
 	itemEnd = decoder->frame.componentList + decoder->frame.componentCount;

    /* Find the maximum factors and the factors for each component. */    
    for (item = decoder->frame.componentList; item < itemEnd; item ++)
    {
        /* Find the opposing scan header component. */
        for (c = 0; ; c ++)
        {
            JPEG_ScanHeader_Component *sc;

            JPEG_Assert (c < decoder->scan.componentCount);
            sc = &decoder->scan.componentList [c];
            if (sc->selector != item->selector)
                continue;
            
            /* Decompress the DC table if necessary. */
            if (sc->dcTable != decoder->dcTableUse[0] && sc->dcTable != decoder->dcTableUse[1])
            {
                const char *tablePointer = decoder->dcTables [(int) sc->dcTable];
                
                if (decoder->dcTableUse [0] == -1)
                    decoder->dcTableUse [0] = sc->dcTable, JPEG_HuffmanTable_Read (&dcTableList [0], &tablePointer);
                else if (decoder->dcTableUse [1] == -1)
                    decoder->dcTableUse [1] = sc->dcTable, JPEG_HuffmanTable_Read (&dcTableList [1], &tablePointer);
                else
                    JPEG_Assert (0);
            }
            
            /* Decompress the AC table if necessary. */
            if (sc->acTable != decoder->acTableUse[0] && sc->acTable != decoder->acTableUse [1])
            {
                const char *tablePointer = decoder->acTables [(int) sc->acTable];
                
                if (decoder->acTableUse [0] == -1)
                    decoder->acTableUse [0] = sc->acTable, JPEG_HuffmanTable_Read (&acTableList [0], &tablePointer);
                else if (decoder->acTableUse [1] == -1)
                    decoder->acTableUse [1] = sc->acTable, JPEG_HuffmanTable_Read (&acTableList [1], &tablePointer);
                else
                    JPEG_Assert (0);
            }
            
            decoder->frameComponents[c] = item;
            break;
        }
        
        /* Add the sum for a later assertion test. */
        decoder->factorSum += item->horzFactor * item->vertFactor;
        
        /* Adjust the maximum horizontal and vertical scaling factors as necessary. */
        if (item->horzFactor > decoder->horzMax)
            decoder->horzMax = item->horzFactor;
        if (item->vertFactor > decoder->vertMax)
            decoder->vertMax = item->vertFactor;
            
        /* Update the relevant component scaling factors if necessary. */
        if (item->selector == 1)
        {
            decoder->YHorzFactor = item->horzFactor;
            decoder->YVertFactor = item->vertFactor;
        }
        else if (item->selector == 2)
        {
            decoder->CbHorzFactor = item->horzFactor;
            decoder->CbVertFactor = item->vertFactor;
        }
        else if (item->selector == 3)
        {
            decoder->CrHorzFactor = item->horzFactor;
            decoder->CrVertFactor = item->vertFactor;
        }
    }
   
    /* Ensure that we have enough memory for these factors. */
    JPEG_Assert (decoder->factorSum < JPEG_MAXIMUM_SCAN_COMPONENT_FACTORS);
     
    /* Split up blockBase according to the components. */
    decoder->YBlock = decoder->blockBase;
    decoder->CbBlock = decoder->YBlock + decoder->YHorzFactor * decoder->YVertFactor * JPEG_DCTSIZE2;
    decoder->CrBlock = decoder->CbBlock + decoder->CbHorzFactor * decoder->CbVertFactor * JPEG_DCTSIZE2;
    
    /* Compute the right shift to be done after multiplying against the scaling factor. */
    if (decoder->horzMax == 1) horzShift = 8;
    else if (decoder->horzMax == 2) horzShift = 7;
    else if (decoder->horzMax == 4) horzShift = 6;
    
    /* Compute the right shift to be done after multiplying against the scaling factor. */
    if (decoder->vertMax == 1) vertShift = 8;
    else if (decoder->vertMax == 2) vertShift = 7;
    else if (decoder->vertMax == 4) vertShift = 6;

    /* Adjust the scaling factors for our parameters. */    
    decoder->YHorzFactor <<= horzShift;
    decoder->YVertFactor <<= vertShift;
    decoder->CbHorzFactor <<= horzShift;
    decoder->CbVertFactor <<= vertShift;
    decoder->CrHorzFactor <<= horzShift;
    decoder->CrVertFactor <<= vertShift;
    
    /* Clear the Cb channel for potential grayscale. */
    {
        signed char *e = decoder->CbBlock + JPEG_DCTSIZE2;
        
        do *-- e = 0;
        while (e > decoder->CbBlock);
    }
    
    /* Clear the Cr channel for potential grayscale. */
    {
        signed char *e = decoder->CrBlock + JPEG_DCTSIZE2;
        
        do *-- e = 0;
        while (e > decoder->CrBlock);
    }

	JPEG_Decoder_ReadImage_Reset();
}

void JPEG_Decoder_ReadImage_Reset()
{
	int c;

	bits_left = bits_data = 0;
    decoder->data = decoder->startdata;
    decoder->restartInterval = decoder->initRestartInterval;
	bytes_left = 1;
	rewind_point = NULL;
    /* Clear the DC parameters. */
    for (c = 0; c < JPEG_MAXIMUM_COMPONENTS; c ++)
        decoder->dcLast[c] = 0;
}

/* Takes information discovered in JPEG_Decoder_ReadHeaders and loads the
 * image.  This is a public function; see gba-jpeg.h for more information on it.
 */
int JPEG_Decoder_ReadImage (JPEG_OUTPUT_TYPE *out, int outWidth, int outHeight)
{
    int c, bx, by, cx, cy; /* Various loop parameters. */
    
	JPEG_Decoder_ReadImage_Reset(decoder);

    int YHorzFactor = decoder->YHorzFactor, YVertFactor = decoder->YVertFactor;
    int CbHorzFactor = decoder->CbHorzFactor, CbVertFactor = decoder->CbVertFactor;
    int CrHorzFactor = decoder->CrHorzFactor, CrVertFactor = decoder->CrVertFactor;
    const char *data = decoder->data;
    signed char *YBlock = decoder->YBlock;
    signed char *CbBlock = decoder->CbBlock;
    signed char *CrBlock = decoder->CrBlock;

	int i, oldline, dxline, line, color;
  
	dxline = oldline = line = 0; 
    /* Now run over each MCU horizontally, then vertically. */
    for (by = 0; by < decoder->frame.height; by += decoder->vertMax * JPEG_DCTSIZE)
    {
		if (by >= outHeight) {
			decoder->data = data;
			return 1;
		}
		// Progress bar
		dxline += 160 * decoder->vertMax * JPEG_DCTSIZE;
		while (dxline >= decoder->frame.height)
		{
			line++;
			dxline -= decoder->frame.height;
		}
		//line = by * 160 / decoder->frame.height;
		for (; oldline < line; oldline++)
		{
		    color = (oldline>>3);
		    color |= (color<<5) | (color<<10);
		    color |= (color<<16);
	    	for (i = 0; i < 120; i++)
			((unsigned int *) (0x06000000))[oldline*120+i] = color; //0x7fff7fff
		}
        for (bx = 0; bx < decoder->frame.width; bx += decoder->horzMax * JPEG_DCTSIZE)
        {
            /* Read the components for the MCU. */
            for (c = 0; c < decoder->scan.componentCount; c ++)
            {
                JPEG_ScanHeader_Component *sc = &decoder->scan.componentList [c];
                JPEG_FrameHeader_Component *fc = decoder->frameComponents[c];
                JPEG_HuffmanTable *dcTable, *acTable;
                JPEG_FIXED_TYPE *quant = decoder->quantTables [(int) fc->quantTable];
                int stride = fc->horzFactor * JPEG_DCTSIZE;
                signed char *chunk = 0;
                
                dcTable = &dcTableList [sc->dcTable == decoder->dcTableUse [1] ? 1 : 0];
                acTable = &acTableList [sc->acTable == decoder->acTableUse [1] ? 1 : 0];
                
                /* Compute the output chunk. */
                if (fc->selector == 1)
                    chunk = YBlock;
                else if (fc->selector == 2)
                    chunk = CbBlock;
                else if (fc->selector == 3)
                    chunk = CrBlock;
                    
                for (cy = 0; cy < fc->vertFactor * JPEG_DCTSIZE; cy += JPEG_DCTSIZE)
                {
                    for (cx = 0; cx < fc->horzFactor * JPEG_DCTSIZE; cx += JPEG_DCTSIZE)
                    {
                        int start = cx + cy * stride;
                        JPEG_FIXED_TYPE zz [JPEG_DCTSIZE2];

                        /* Decode coefficients. */
                        DecodeCoefficients (&decoder->dcLast[c], zz, quant, dcTable, acTable, &data, ToZigZag);

                        /* Perform an IDCT if this component will contribute to the image. */
                        if (chunk)
                        {
                            IDCT_Columns (zz);
                            IDCT_Rows (zz, chunk + start, stride);
                        }
                    }
                }
            }
			if (bx >= outWidth)
				continue;
            /* Check that our block will be in-range; this should actually use clamping. */
            // if (bx + horzMax * JPEG_DCTSIZE > outWidth || by + vertMax * JPEG_DCTSIZE > outHeight)
            //    continue;
                
            /* Convert our block from YCbCr to the output. */
            ConvertBlock (YBlock, CbBlock, CrBlock,
                YHorzFactor, YVertFactor, CbHorzFactor, CbVertFactor, CrHorzFactor, CrVertFactor,
                decoder->horzMax * JPEG_DCTSIZE, decoder->vertMax * JPEG_DCTSIZE, out + bx + by * outWidth, bx, by, outWidth, outHeight, ComponentRange);
            
            /* Handle the restart interval. */
            if (decoder->initRestartInterval && --decoder->restartInterval == 0)
            {
                decoder->restartInterval = decoder->initRestartInterval;
                JPEG_BITS_REWIND ();
                if (((data [0] << 8) | data [1]) == JPEG_Marker_EOI)
                    goto finish;
                JPEG_Assert (data [0] == 0xFF && (data [1] >= 0xD0 && data [1] <= 0xD7));
                for (c = 0; c < JPEG_MAXIMUM_COMPONENTS; c ++)
                    decoder->dcLast[c] = 0;
                data += 2;
            }
        }
    }
   
finish:
    /* Make sure we read an EOI marker. */ 
    JPEG_BITS_REWIND ();
    JPEG_Assert (((data [0] << 8) | data [1]) == JPEG_Marker_EOI);
    data += 2;
    
    /* Clear up and return success. */
    decoder->data = data;
    return 1;
}

/* Read an JPEG_Marker_SOFn marker into frame.  This expects to start
 * processing immediately after the marker.
 */
int JPEG_FrameHeader_Read (JPEG_FrameHeader *frame, const char **dataBase, JPEG_Marker marker)
{
    const char *data = *dataBase;
    unsigned short length = (data [0] << 8) | data [1];
    int index;

    (void) length;
    JPEG_Assert (length >= 8);        
    data += 2; /* Skip the length. */
    frame->marker = marker;
    frame->encoding = (marker >= 0xFFC0 && marker <= 0xFFC7) ? 0 : 1;
    frame->differential = !(marker >= 0xFFC0 && marker <= 0xFFC3 && marker >= 0xFFC8 && marker <= 0xFFCB);
    
    frame->precision = *data ++;
    frame->height = (data [0] << 8) | data [1]; data += 2;
    frame->width = (data [0] << 8) | data [1]; data += 2;
    frame->componentCount = *data ++;

	// leave a permanent record of picture size
	jpegHeight = frame->height;
	jpegWidth = frame->width;
    

    JPEG_Assert (frame->precision == 8);
    JPEG_Assert (frame->componentCount <= JPEG_MAXIMUM_COMPONENTS);
    JPEG_Assert (length == 8 + 3 * frame->componentCount);
    
    /* Read the frame components. */
    for (index = 0; index < frame->componentCount; index ++)
    {
        JPEG_FrameHeader_Component *c = &frame->componentList [index];
        char pair;
        
        c->selector = *data ++;
        pair = *data ++;
		if (frame->componentCount == 1) {
			c->horzFactor = c->vertFactor = 1;
		} else {
	        c->horzFactor = pair >> 4;
    	    c->vertFactor = pair & 15;
		}
        c->quantTable = *data ++;
        
        JPEG_Assert (c->horzFactor == 1 || c->horzFactor == 2 || c->horzFactor == 4);
        JPEG_Assert (c->vertFactor == 1 || c->vertFactor == 2 || c->vertFactor == 4);
        JPEG_Assert (c->quantTable <= 3);
    }
    
    *dataBase = data;
    return 1;
}

/* Read a JPEG_Marker_SOS marker into scan.  This expects to start processing
 * immediately after the marker.
 */
int JPEG_ScanHeader_Read (JPEG_ScanHeader *scan, const char **dataBase)
{
    const char *data = *dataBase;
    unsigned short length = (data [0] << 8) | data [1];
    JPEG_ScanHeader_Component *c, *cEnd;
    char pair;
 
    (void) length;
    JPEG_Assert (length >= 6);
    data += 2; /* Skip the length. */
    scan->componentCount = *data ++;
    
    JPEG_Assert (scan->componentCount <= JPEG_MAXIMUM_COMPONENTS);
    JPEG_Assert (length == 6 + 2 * scan->componentCount);

    /* Read the scan components. */    
    for (c = scan->componentList, cEnd = c + scan->componentCount; c < cEnd; c ++)
    {
        c->selector = *data ++;
        pair = *data ++;
        c->dcTable = pair >> 4;
        c->acTable = pair & 15;
        
        JPEG_Assert (c->dcTable < 4);
        JPEG_Assert (c->acTable < 4);
    }
    
    /* Read the spectral and approximation footers, which are used for
     * progressive.
     */
     
    scan->spectralStart = *data ++;
    scan->spectralEnd = *data ++;
    JPEG_Assert (scan->spectralStart <= 63);
    JPEG_Assert (scan->spectralEnd <= 63);
    pair = *data ++;
    scan->successiveApproximationBitPositionHigh = pair >> 4;
    scan->successiveApproximationBitPositionLow = pair & 15;
    JPEG_Assert (scan->successiveApproximationBitPositionHigh <= 13);
    JPEG_Assert (scan->successiveApproximationBitPositionLow <= 15);
    
    *dataBase = data;
    return 1;
}

/* Read all headers from the very start of the JFIF stream to right after the
 * SOS marker.
 */
 
int JPEG_Decoder_ReadHeaders (JPEG_Decoder *decoder, const char **dataBase)
{
    const char *data = *dataBase;
    JPEG_Marker marker;
    int c;
 
    /* Initialize state and assure that this is a JFIF file. */   
    decoder->initRestartInterval = 0;
    JPEG_Assert (((data [0] << 8) | data [1]) == JPEG_Marker_SOI);
    data += 2;
    
    /* Start reading every marker as it comes in. */
    while (1)
    {
        marker = (JPEG_Marker)((data [0] << 8) | data [1]);
        data += 2;
        
        switch (marker)
        {
            /* This block is just skipped over. */
            case JPEG_Marker_APP0:
            case JPEG_Marker_APP1:
            case JPEG_Marker_APP2:
            case JPEG_Marker_APP3:
            case JPEG_Marker_APP4:
            case JPEG_Marker_APP5:
            case JPEG_Marker_APP6:
            case JPEG_Marker_APP7:
            case JPEG_Marker_APP8:
            case JPEG_Marker_APP9:
            case JPEG_Marker_APP10:
            case JPEG_Marker_APP11:
            case JPEG_Marker_APP12:
            case JPEG_Marker_APP13:
            case JPEG_Marker_APP14:
            case JPEG_Marker_APP15:
            case JPEG_Marker_COM:
                data += (data [0] << 8) | data [1];
                break;
            
            case JPEG_Marker_DHT: /* Define Huffman table.  We just skip it for later decompression. */
            {
                unsigned short length = (data [0] << 8) | data [1];
                const char *end = data + length;
                
                JPEG_Assert (length >= 2);
                data += 2;
                while (data < end)
                {
                    char pair, type, slot;
                    
                    pair = *data ++;
                    type = pair >> 4;
                    slot = pair & 15;
                    
                    JPEG_Assert (type == 0 || type == 1);
                    JPEG_Assert (slot <= 15);
                    
                    if (type == 0)
                        decoder->dcTables [(int) slot] = data;
                    else
                        decoder->acTables [(int) slot] = data;
                        
                    if (!JPEG_HuffmanTable_Skip (&data))
                        return 0;
                }
                
                JPEG_Assert (data == end);
                break;
            }
            
            case JPEG_Marker_DQT: /* Define quantization table. */
            {
                unsigned short length = (data [0] << 8) | data [1];
                const char *end = data + length;
                int col, row;
                JPEG_FIXED_TYPE *s;
                
                JPEG_Assert (length >= 2);
                data += 2;
                
                while (data < end)
                {
                    int pair, slot, precision;
                    
                    pair = *data ++;
                    precision = pair >> 4;
                    slot = pair & 15;
                    
                    JPEG_Assert (precision == 0); /* Only allow 8-bit. */
                    JPEG_Assert (slot < 4); /* Ensure the slot is in-range. */
                    JPEG_Assert (data + 64 <= end); /* Ensure it's the right size. */
                    
                    s = decoder->quantTables [slot];
                   
                    for (c = 0; c < JPEG_DCTSIZE2; c ++)
                        s [c] = JPEG_ITOFIX (*data ++);
                    
                    /* Multiply against the AAN factors. */
                    for (row = 0; row < JPEG_DCTSIZE; row ++)
                        for (col = 0; col < JPEG_DCTSIZE; col ++)
                        {
                            JPEG_FIXED_TYPE *item = &s [col + row * JPEG_DCTSIZE];
                            
                            *item = JPEG_FIXMUL (*item, JPEG_AANScaleFactor [(int) ToZigZag [row * JPEG_DCTSIZE + col]]);
                        }
                }
                
                JPEG_Assert (data == end); /* Ensure we've finished it. */
                break;
            }
        
            case JPEG_Marker_DRI: /* Define restart interval. */
                JPEG_Assert (((data [0] << 8) | data [1]) == 4); /* Check the length. */
                decoder->initRestartInterval = (data [2] << 8) | data [3];
                data += 4;
                break;
            
            case JPEG_Marker_SOF0: /* Start of Frame: Baseline Sequential Huffman. */
                if (!JPEG_FrameHeader_Read (&decoder->frame, &data, marker))
                    return 0;
                break;
            
            case JPEG_Marker_SOS: /* Start of scan, immediately followed by the image. */
                if (!JPEG_ScanHeader_Read (&decoder->scan, &data))
                    return 0;
                *dataBase = data;
                return 1;
                
            default: /* No known marker of this type. */
                JPEG_Assert (0);
                break;
        }
    }
}

/* Skip past a Huffman table section.  This expects to be called after reading
 * the DHT marker and the type/slot pair.
 */
int JPEG_HuffmanTable_Skip (const char **dataBase)
{
    const char *data = *dataBase;
    int c, total = 16;
    
    for (c = 0; c < 16; c ++)
        total += *data ++;
    *dataBase += total;
    return 1;
}

/* Decode a Huffman table and initialize its data.  This expects to be called
 * after the DHT marker and the type/slot pair.
 */
int JPEG_HuffmanTable_Read (JPEG_HuffmanTable *huffmanTable, const char **dataBase)
{
    const char *data = *dataBase;
    const char *bits;
    int huffcode [256];
    char huffsize [256];
    int total = 0;
    int c;
    
    bits = data;
    for (c = 0; c < 16; c ++)
        total += *data ++;
    huffmanTable->huffval = data;
    data += total;
    
    /*void GenerateSizeTable ()*/
    {
        int k = 0, i = 1, j = 1;
        
        do
        {
            while (j ++ <= bits [i - 1])
                huffsize [k ++] = i;
            i ++;
            j = 1;
        }
        while (i <= 16);
            
        huffsize [k] = 0;
    }
    
    /*void GenerateCodeTable ()*/
    {
        int k = 0, code = 0, si = huffsize [0];

        while (1)
        {            
            do huffcode [k ++] = code ++;
            while (huffsize [k] == si);
                
            if (huffsize [k] == 0)
                break;
            
            do code <<= 1, si ++;
            while (huffsize [k] != si);
        }
    }
    
	// Ensure the huffman decoder terminates
	huffmanTable->maxcode[16] = 0xfffff;
    /*void DecoderTables ()*/
    {
        int i = 0, j = 0;
        
        while (1)
        {
            if (i >= 16)
                break;
            if (bits [i] == 0)
                huffmanTable->maxcode [i] = -1;
            else
            {
                huffmanTable->valptr [i] = &huffmanTable->huffval [j - huffcode [j]];
                j += bits [i];
                huffmanTable->maxcode [i] = huffcode [j - 1];
            }
            i ++;
        }
    }
    
    /*void GenerateLookahead ()*/
    {
        int l, i, p, c, ctr;
        
        for (c = 0; c < 256; c ++)
            huffmanTable->look_nbits [c] = 0;
            
        p = 0;
        for (l = 1; l <= 8; l ++)
        {
            for (i = 1; i <= bits [l - 1]; i ++, p ++)
            {
                int lookbits = huffcode [p] << (8 - l);
                
                for (ctr = 1 << (8 - l); ctr > 0; ctr --)
                {
                    huffmanTable->look_nbits [lookbits] = l;
                    huffmanTable->look_sym [lookbits] = huffmanTable->huffval [p];
                    lookbits ++;
                }
            }
        }
    }
    
    *dataBase = data;
    return 1;
}

#define MIN(a,b) (a<b ? a : b)

const char *jpg_data;

/* Perform the two steps necessary to decompress a JPEG image.
 * Nothing fancy about it.
 */
int JPEG_DecompressImage_Init (JPEG_Decoder *dec, const char *data, JPEG_OUTPUT_TYPE **out, int *outWidth, int *outHeight)
{
	int space_left;

	decoder = dec;
	acTableList = pmalloc(sizeof(JPEG_HuffmanTable)*2);
	dcTableList = pmalloc(sizeof(JPEG_HuffmanTable)*2);

	space_left = pmemory_free();
	//fprintf(stderr, "space left: %d\n", space_left);

    // Clear memory.
    /*for (i = 0; i < 64*1024; i++)
	    ((int *)(0x02000000))[i] = 0; */
	if (!JPEG_Match(data))
		return 0;
    
    if (!JPEG_Decoder_ReadHeaders (decoder, &data))
        return 0;

    if (jpegWidth * jpegHeight * 2 > space_left)
		*outHeight = space_left/jpegWidth/2;
    else
	    *outHeight = jpegHeight;
    
	*outWidth = jpegWidth;

    *out = pmalloc(*outWidth * *outHeight * 2);
    if (!*out)
	    return 2;

	jpg_data = data;

	JPEG_Decoder_ReadImage_Init(jpg_data);

    return 1;
}

int JPEG_DecompressImage(JPEG_OUTPUT_TYPE *out, int outWidth, int outHeight)
{
	if (!JPEG_Decoder_ReadImage(out, outWidth, outHeight))
		return 0;
	return 1;
}

/* Return whether this code is a JPEG file.  Unfortunately it will incorrectly
 * match variants such as JPEG 2000 and JPEG-LS.  A better function would
 * skip known markers until it reaches an unknown marker or a handled
 * SOFn.
 */
 
int JPEG_Match (const char *data)
{
    if (data [0] != 0xFF) return 0;
    if (data [1] != 0xD8) return 0;
    if (data [2] != 0xFF) return 0;
    if (data [3] != 0xE0) return 0;
	if (data [6] != 0x4A) return 0;
	if (data [7] != 0x46) return 0;
	if (data [8] != 0x49) return 0;
	if (data [9] == 0x46) return 1;
    return 0;
}
