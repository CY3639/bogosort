/*
 * Implementation of the shared interface declared in common.h.
 *
 * This file is compiled into every executable:
 *   gcc -O2 -o greedy  common.c greedy.c
 *   gcc -O2 -o exact   common.c exact.c
 *
 * No main() lives here. Only infrastructure that both solvers share.
 */

#include "common.h"

#include <stdio.h>     // getline(), fopen(), fclose(), fprintf(), fputs(), fputc() 
#include <stdlib.h>    // malloc(), calloc(), free()                                
#include <string.h>    // strcmp(), strncmp(), strstr(), strlen()                   
#include <assert.h>    // assert()                                                  


/* 
 * Section 2 - Fragment I/O
 */

Fragment *
read_all_fragments( char const * file_name )
{
    FILE * input = stdin;
    if( strcmp( file_name, "-" ) != 0 )
        input = fopen( file_name, "r" );

    if( input == NULL ) {
        fprintf( stderr, "Error: file could not be opened for reading: %s\n", file_name );
        return NULL;
    }

    Fragment * top_fragment = NULL;

    while( 1 ) {
        char   * line_buffer = NULL;
        size_t   buffer_size = 0;

        ssize_t read_bytes = getline( &line_buffer, &buffer_size, input );
        if( read_bytes <= 0 ) {
            if( line_buffer != NULL ) free( line_buffer );
            break;
        }

        assert( line_buffer != NULL );
        assert( line_buffer[read_bytes] == '\0' );

        /* Strip trailing newline / carriage-return. */
        for( int i = 0; i < read_bytes; i++ )
            if( line_buffer[i] == '\n' || line_buffer[i] == '\r' )
                line_buffer[i] = '\0';

        /* Skip blank lines — they carry no fragment information. */
        if( line_buffer[0] == '\0' ) {
            free( line_buffer );
            continue;
        }

        Fragment * new_fragment = malloc( sizeof( Fragment ) );
        assert( new_fragment != NULL );
        new_fragment->next_fragment  = top_fragment;
        new_fragment->fragment_string = line_buffer;
        line_buffer   = NULL;   /* ownership transferred */
        top_fragment  = new_fragment;
    }

    fclose( input );
    return top_fragment;
}

void
free_all_fragments( Fragment * top_fragment )
{
    while( top_fragment != NULL ) {
        Fragment * this_fragment = top_fragment;
        top_fragment = this_fragment->next_fragment;
        free( this_fragment->fragment_string );
        free( this_fragment );
    }
}


/* 
 * Section 3 - Preprocessing
 */

Fragment *
remove_substring_fragments( Fragment * top_fragment )
{
    Fragment * prev    = NULL;
    Fragment * current = top_fragment;

    while( current != NULL ) {
        int removed = 0;

        // Check if any other fragment contains current as a substring.
        for( Fragment * other = top_fragment; other != NULL; other = other->next_fragment ) {
            if( other == current ) continue;
            if( strstr( other->fragment_string, current->fragment_string ) != NULL ) {
                /* current is redundant — remove it from the list. */
                Fragment * to_free = current;
                if( prev == NULL ) {
                    top_fragment = current->next_fragment;
                    current      = top_fragment;
                } else {
                    prev->next_fragment = current->next_fragment;
                    current             = current->next_fragment;
                }
                free( to_free->fragment_string );
                free( to_free );
                removed = 1;
                break;
            }
        }

        if( !removed ) {
            prev    = current;
            current = current->next_fragment;
        }
    }

    return top_fragment;
}

int **
build_overlap_matrix( Fragment * top_fragment, int * fragment_count )
{
    /* Count fragments. */
    int n = 0;
    for( Fragment * p = top_fragment; p != NULL; p = p->next_fragment ) n++;

    /* Allocate n x n matrix, zero-initialised. */
    int ** matrix = malloc( n * sizeof( int * ) );
    assert( matrix != NULL );
    for( int i = 0; i < n; i++ ) {
        matrix[i] = calloc( n, sizeof( int ) );
        assert( matrix[i] != NULL );
    }

    /*
     * For each ordered pair (i, j), find the longest k such that
     * the last k characters of fragment[i] equal the first k characters
     * of fragment[j].  That is the overlap weight for edge i -> j.
     */
    int i = 0;
    for( Fragment * fi = top_fragment; fi != NULL; fi = fi->next_fragment, i++ ) {
        size_t len_i = strlen( fi->fragment_string );
        int j = 0;
        for( Fragment * fj = top_fragment; fj != NULL; fj = fj->next_fragment, j++ ) {
            if( fi == fj ) continue;
            size_t len_j   = strlen( fj->fragment_string );
            size_t max_ov  = len_i < len_j ? len_i : len_j;

            for( size_t k = max_ov; k > 0; k-- ) {
                /* suffix of fi (length k) == prefix of fj (length k)? */
                if( strncmp( fi->fragment_string + len_i - k,
                             fj->fragment_string,
                             k ) == 0 ) {
                    matrix[i][j] = (int)k;
                    break;
                }
            }
        }
    }

    *fragment_count = n;
    return matrix;
}


/* 
 * Section 4 - Fragment array utility
 */

Fragment **
fragments_to_array( Fragment * top_fragment, int count )
{
    Fragment ** arr = malloc( count * sizeof( Fragment * ) );
    assert( arr != NULL );
    Fragment * p = top_fragment;
    for( int i = 0; i < count; i++, p = p->next_fragment )
        arr[i] = p;
    return arr;
}


/* 
 * Section 5 - Solution validation
 */

bool
solution_is_valid( char const * candidate, Fragment * top_fragment )
{
    for( Fragment * f = top_fragment; f != NULL; f = f->next_fragment ) {
        if( strstr( candidate, f->fragment_string ) == NULL )
            return false;
    }
    return true;
}


/* 
 * Section 6 - Output helpers
 */

void
print_solution( char const * solution )
{
    fputs( solution, stdout );
    fputc( '\n',    stdout );
}