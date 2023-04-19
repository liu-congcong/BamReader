#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <limits.h>
#include <zlib.h>
#include <inttypes.h>

#define BUFFER 1024 * 1024

typedef struct READ
{
    int32_t position;
    uint32_t *cigar;
    uint8_t *read;
    struct READ *next;
} READ;

typedef struct
{
    char *sequence_name;
    int32_t *sequence_length;
    READ *read;
} SEQUENCE;

int free_sequence_list(SEQUENCE *sequence_list, int32_t sequences)
{
    for (int32_t sequence_index = 0; sequence_index < sequences; sequence_index++)
    {
        free(sequence_list[sequence_index].sequence_name);
    }
    free(sequence_list);
    return 0;
}

bool is_sorted_file(char *header)
{
    *strchr(header, '\n') = 0;
    if (!strncmp(header, "@HD", 3) && strstr(header, "SO:coordinate"))
        return true;
    else
        return false;
}

int read_header(gzFile open_file, SEQUENCE **sequence_list, int32_t *sequences)
{
    char magic[4];
    gzread(open_file, magic, sizeof(char) * 4);
    assert(!strncmp(magic, "BAM\1", 4)); /* This is a BAM file. */
    int32_t header_length;
    gzread(open_file, &header_length, sizeof(int32_t));
    char *header = malloc(sizeof(char) * header_length);
    assert(header);
    gzread(open_file, header, sizeof(char) * header_length); /* @ lines. */
    /*
    if (!is_sorted_file(header))
    {
        printf("BAM file needs be sorted.\n");
        exit(EXIT_FAILURE);
    }
    */
    free(header); /* ... */
    gzread(open_file, sequences, sizeof(int32_t));
    *sequence_list = malloc(sizeof(SEQUENCE) * (*sequences));
    assert(*sequence_list);
    for (int32_t sequence_index = 0; sequence_index < (*sequences); sequence_index++) /* Read sequence list. */
    {
        int32_t sequence_name_length;
        gzread(open_file, &sequence_name_length, sizeof(int32_t));
        (*sequence_list)[sequence_index].sequence_name = malloc(sizeof(char) * sequence_name_length);
        gzread(open_file, (*sequence_list)[sequence_index].sequence_name, sizeof(char) * sequence_name_length);
        int32_t sequence_length;
        gzread(open_file, &(*sequence_list)[sequence_index].sequence_length, sizeof(int32_t));
    }
    return 0;
}

int decode_sequence(uint8_t *sequence, int32_t sequence_length, char *decoded_sequence)
{
    char decoder[16] = "=ACMGRSVTWYHKDBN";
    decoded_sequence[sequence_length * 2 + 1] = 0;
    for (int32_t index = 0; index < sequence_length; index++)
    {
        decoded_sequence[index * 2] = decoder[sequence[index] >> 4];
        decoded_sequence[index * 2 + 1] = decoder[sequence[index] & 0x0F];
    }
    return 0;
}

int decode_cigar_debug(uint32_t *cigar, uint16_t operations, char *decoded_cigar, char *buffer)
{
    char decoder[9] = "MIDNSHP=X";
    buffer[0] = 0;
    decoded_cigar[0] = 0;
    for (uint16_t operation = 0; operation < operations; operation++)
    {
        sprintf(buffer, "%d%c", cigar[operation] >> 4, decoder[cigar[operation] & 0x0000000f]);
        strcat(decoded_cigar, buffer);
    }
    decoded_cigar[strlen(decoded_cigar)] = 0;
    return 0;
}

int add_alignment()
{
    return 0;
}

int read_alignments(gzFile open_file, SEQUENCE *sequence_list, uint8_t min_mapq)
{
    char *buffer1 = malloc(sizeof(char) * BUFFER);
    assert(buffer1);
    uint32_t *buffer2 = malloc(sizeof(uint32_t) * UINT16_MAX);
    assert(buffer2);
    uint8_t *buffer3 = malloc(sizeof(uint8_t) * BUFFER);
    assert(buffer3);

    char *buffer4 = malloc(sizeof(char) * BUFFER);
    assert(buffer4);
    char *buffer5 = malloc(sizeof(char) * BUFFER);
    assert(buffer5);

    int32_t sequences = 0;
    int32_t alignment_length = 0;
    int32_t sequence_index;
    int32_t position;
    uint8_t read_name_length;
    uint8_t mapq;
    uint16_t bai_bin;
    uint16_t cigar_operations;
    uint16_t flag;
    int32_t read_length;
    int32_t next_sequence_index;
    int32_t next_position;
    int32_t template_length;
    printf("Sequence ID\tPosition\tCigar\tRead\n"); /* header */
    while (gzread(open_file, &alignment_length, sizeof(int32_t))) /* Return size_t. */
    {
        gzread(open_file, &sequence_index, sizeof(int32_t)); /* refID */
        gzread(open_file, &position, sizeof(int32_t)); /* pos */
        gzread(open_file, &read_name_length, sizeof(uint8_t)); /* l_read_name */
        gzread(open_file, &mapq, sizeof(uint8_t)); /* mapq */
        gzread(open_file, &bai_bin, sizeof(uint16_t)); /* bin */
        gzread(open_file, &cigar_operations, sizeof(uint16_t)); /* n_cigar_op */
        gzread(open_file, &flag, sizeof(uint16_t)); /* flag */
        gzread(open_file, &read_length, sizeof(int32_t)); /* l_seq */
        gzread(open_file, &next_sequence_index, sizeof(int32_t)); /* next_refID */
        gzread(open_file, &next_position, sizeof(int32_t)); /* next_pos */
        gzread(open_file, &template_length, sizeof(int32_t)); /* tlen */
        gzread(open_file, buffer1, sizeof(char) * read_name_length); /* read_name */
        gzread(open_file, buffer2, sizeof(uint32_t) * cigar_operations); /* cigar */
        gzread(open_file, buffer3, sizeof(uint8_t) * (read_length + 1) / 2); /* seq: =ACMGRSVTWYHKDBN: A:1, C:2, G:4, T:8 */
        gzread(open_file, buffer1, sizeof(char) * read_length); /* qual */
        alignment_length -= sizeof(int32_t) * 6 + sizeof(uint8_t) * (2 + (read_length + 1) / 2) + sizeof(uint16_t) * 3 + sizeof(char) * (read_name_length + read_length) + sizeof(uint32_t) * cigar_operations;
        gzread(open_file, buffer1, sizeof(char) * alignment_length); /* auxiliary data */
        if (sequence_index >= 0 && mapq >= min_mapq)
        {
            decode_cigar_debug(buffer2, cigar_operations, buffer4, buffer1);
            decode_sequence(buffer3, (read_length + 1) / 2, buffer5);
            // add_alignment(sequence_list, sequence_index, position);
            printf("%s\t%d\t%s\t%s\n",sequence_list[sequence_index].sequence_name, position + 1, buffer4, buffer5);
        }
    }
    free(buffer1);
    free(buffer2);
    free(buffer3);
    free(buffer4);
    free(buffer5);
    return 0;
}

int read_bam(char *bam, uint8_t min_mapq)
{
    SEQUENCE *sequence_list = NULL;
    int32_t sequences = 0;
    gzFile open_file = gzopen(bam, "rb");
    read_header(open_file, &sequence_list, &sequences);
    /*
    for (int32_t sequence_index = 0; sequence_index < sequences; sequence_index++)
    {
        printf("sequence id: %s, sequence length: %d\n", sequence_list[sequence_index].sequence_name, sequence_list[sequence_index].sequence_length);
    }
    */
    assert(sequences); /* Have no sequences. */
    read_alignments(open_file, sequence_list, min_mapq);
    gzclose(open_file);
    free_sequence_list(sequence_list, sequences);
    return 0;
}

int main(int argc, char *argv[])
{
    char input_bam[FILENAME_MAX];
    uint8_t MAPQ = 20;
    int necessary_parameters = 0;
    for (int index = 0; index < argc; index++)
    {
        if (!strncasecmp(argv[index], "-b", 2) || !strncasecmp(argv[index], "--b", 3))
        {
            strncpy(input_bam, argv[index + 1], FILENAME_MAX);
            input_bam[FILENAME_MAX - 1] = 0;
            necessary_parameters++;
        }
        else if (!strncasecmp(argv[index], "-m", 2) || !strncasecmp(argv[index], "--m", 3))
        {
            sscanf(argv[index + 1], "%"SCNu8, &MAPQ);
        }
    }
    if (necessary_parameters == 1)
        read_bam(input_bam, MAPQ);
    else
    {
        printf("Usage: %s -bam SORTED.BAM -mapq MAPQ > OUTPUT\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    return 0;
}
