/*
 * Shared interface for Minimum-Size Text Reconstruction from Overlapping Fragments.
 *
 * Example: Neither greedy.c nor exact.c should duplicate the declarations below.
 * Compile each executable against common.c, e.g.:
 *   gcc -O2 -o greedy  common.c greedy.c
 *   gcc -O2 -o exact   common.c exact.c
 * 
 * this satisfies the project specification for "one or several self-contained C programs" requirement,
 * while avoiding duplicated logic.
 *
 * Any change to this file must be agreed upon by all group members,
 * since it affects every executable simultaneously.
 * 
 */

#ifndef COMMON_H
#define COMMON_H

#include <stddef.h>    // size_t                
#include <stdbool.h>   // bool                    


/* 
 * Section 1 - Core data structure
 *
 * A singly linked list node representing one input fragment.
 * callee frees via free_all_fragments().
 */

typedef struct fragment_s {
    struct fragment_s * next_fragment;
    char              * fragment_string;   // null-terminated; heap-allocated
} Fragment;


/* 
 * Section 2 - Fragment I/O
 *
 * Reads every non-empty line from `file_name` (or stdin if "-") and returns
 * them as a singly linked list of Fragment nodes.
 *
 * caller must eventually call free_all_fragments() on the result.
 * Returns NULL and prints to stderr on I/O error.
 *
 * Complexity: O(L)  where L = total characters across all fragments.
 */
Fragment * read_all_fragments( char const * file_name );

/*
 * Releases every node and its fragment_string in the list.
 * Safe to call with NULL.
 *
 * Complexity: O(n)  where n = number of fragments.
 */
void free_all_fragments( Fragment * top_fragment );


/* 
 * Section 3 - Preprocessing
 *
 * These functions reduce the problem size before any search begins.
 * Both the greedy and exact solvers should call them in this order.
 */

/*
 * Removes any fragment that is already a substring of another fragment in
 * the list. Such fragments are redundant: any superstring containing the
 * longer fragment automatically contains the shorter one.
 *
 * Returns the (possibly shorter) list head; the original pointer may be
 * invalidated if the first node is removed.
 *
 * Complexity: O(n^2 * m^2)
 *   n = number of fragments, m = length of longest fragment.
 */
Fragment * remove_substring_fragments( Fragment * top_fragment );

/*
 * Computes, for every ordered pair (i, j) of remaining fragments, the length
 * of the longest suffix of fragment[i] that is a prefix of fragment[j].
 * This is the "overlap" weight for the directed edge i -> j in the overlap
 * graph.
 *
 * Writes the fragment count into *fragment_count.
 * Returns a heap-allocated n x n matrix; caller must free each row then the
 * outer pointer.
 *
 * Complexity: O(n^2 * m^2)
 */
int ** build_overlap_matrix( Fragment * top_fragment, int * fragment_count );


/* 
 * Section 4 - Fragment array utility
 *
 * Both solvers need random access to fragments by index. 
 * This converts the linked list to a plain array for O(1) index lookup.
 * 
 * Returns a heap-allocated array of Fragment pointers in list order.
 * `count` must equal the length of the list (i.e. the value written by
 * build_overlap_matrix into fragment_count).
 *
 * Caller frees the array itself (not the nodes, which remain in the list).
 *
 * Complexity: O(n)
 */
Fragment ** fragments_to_array( Fragment * top_fragment, int count );


/* 
 * Section 5 - Solution validation
 *
 * A correct reconstruction must contain every input fragment as a substring.
 * Both solvers should use this to sanity-check output before printing.
 * 
 * Returns true iff every fragment in the list appears as a substring of
 * `candidate`.
 *
 * Complexity: O(n * |candidate| * m)
 */
bool solution_is_valid( char const * candidate, Fragment * top_fragment );


/* 
 * Section 6 - Output helpers
 *
 * Writes `solution` followed by a newline to stdout.
 * Centralised here so both executables format output identically.
 */
void print_solution( char const * solution );


#endif /* COMMON_H */