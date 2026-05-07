/* compile with:
    gcc -O2 -o mgreedy_reconstruction common.c mgreedy.c
*/

/* Overall complexity:
 * MGREEDY phases:   O(n^2 m^2) for overlap matrix (done in common.c)
 *                   O(n^2 log n) for edge sorting + greedy selection
 * 
 * Held-Karp algorithm
 * DP if n < 18:     O(2^n * n^2) time, O(2^n * n) space
 * n can be changed  2^18 = 262,144; table entries = 4,718,592; each entry = 4 bytes
 * considering       4718592*4 = ~18MB RAM
 * available RAM     2^22 = 4,194,304; table entries = 92,274,688; ~360MB
 *                   2^25 = 33,554,432; table entries = 838,860,800; ~3.2GB
*/

#include "common.h"
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

// Section 1 - Edge structure and helpers

typedef struct {
    int from;
    int to;
    int weight; // overlap length
} Edge;
 
static int
edge_cmp_desc(const void *a, const void *b)
{
    return ((const Edge *)b)->weight - ((const Edge *)a)->weight;
}

/* Section 2 - MGREEDY algorithm for optimal cycle cover using greedy edge picking
*  
*  Identical to the greedy algorithm used in kaito.c except edges that close cycles are accepted.
*  Only enforce:
*  - each vertex has at most 1 outgoing selected edge
*  - each vertex has at most 1 incoming selected edge
*  to produce the maximum weight cycle cover of the overlap graph
*/ 

/*
 * Returns an n*n matrix where selected[i][j] = overlap weight
 * if edge i -> j was selected, 0 otherwise.
 *
 * Complexity: O(n^2 log n) for sorting + scanning edges.
 */

static int **
mgreedy_cycle_cover(int **overlap, int n)
{
    // collect all positive weight edges
    int capacity = n * n;
    Edge *edges = malloc(capacity * sizeof(Edge));
    assert(edges != NULL);
    int edge_count = 0;
 
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j && overlap[i][j] > 0) {
                edges[edge_count].from = i;
                edges[edge_count].to = j;
                edges[edge_count].weight = overlap[i][j];
                edge_count++;
            }
        }
    }
 
    qsort(edges, edge_count, sizeof(Edge), edge_cmp_desc);
 
    // track per-vertex in-degree and out-degree (max 1 each)
    int *has_out = calloc(n, sizeof(int)); // has_out[i] = 1 if i already has an outgoing edge
    int *has_in  = calloc(n, sizeof(int)); // has_in[j]  = 1 if j already has an incoming edge
    assert(has_out != NULL && has_in != NULL);
 
    int **selected = malloc(n * sizeof(int *));
    assert(selected != NULL);
    for (int i = 0; i < n; i++) {
        selected[i] = calloc(n, sizeof(int));
        assert(selected[i] != NULL);
    }
 
    for (int e = 0; e < edge_count; e++) {
        int u = edges[e].from;
        int v = edges[e].to;
 
        // degree constraint only, no cycle check.
        if (has_out[u] || has_in[v])
            continue;
 
        selected[u][v] = edges[e].weight;
        has_out[u] = 1;
        has_in[v] = 1;
    }
 
    free(edges);
    free(has_out);
    free(has_in);
    return selected;
}

/* Section 3 - open cycles to remove weakest edge and build chains
*  
*  each cycle in the cycle cover has a cycle-closing edge
*  the one with smallest overlap, the last one MGREEDY added
*  removing it turns the cycle into a linear chain
*  then concatenate each chain into a representative string
*/ 

/*
 * Finds the successor of vertex `v` in the selected edge matrix.
 * Returns -1 if v has no outgoing selected edge.
 */

static int
successor(int **selected, int n, int v)
{
    for (int j = 0; j < n; j++) {
        if (selected[v][j] > 0)
            return j;
    }
    return -1;
}

/* build representative strings from the cycle cover
*
*  - identify all cycles by following successor links
*  - for each cycle, find the edge with minimum overlap (cycle-closing)
*  - remove that edge, the cycle becomes a chain
*  - walk to the chain, merging fragments using their overlaps
*
*  returns a heap allocated array of representative strings
*/

static char **
open_cycles_to_representatives(int **selected, int n, Fragment **frag_arr, int *rep_count)
{
    int *visited = calloc(n, sizeof(int));
    assert(visited != NULL);
 
    char **reps = malloc(n * sizeof(char *)); // at most n representatives
    assert(reps != NULL);
    int count = 0;
 
    for (int start = 0; start < n; start++) {
        if (visited[start])
            continue;
 
        // Trace the cycle starting from start
        int cycle_len = 0;
        int *cycle_nodes = malloc(n * sizeof(int));
        assert(cycle_nodes != NULL);
 
        int cur = start;
        while (!visited[cur]) {
            visited[cur] = 1;
            cycle_nodes[cycle_len++] = cur;
            int next = successor(selected, n, cur);
            if (next == -1)
                break; // isolated node or end of a chain, shouldn't happen in a cycle cover but safeguard
            cur = next;
        }
 
        if (cycle_len == 1 && successor(selected, n, start) == -1) {
            // isolated node - no outgoing or incoming edges selected
            reps[count] = malloc(strlen(frag_arr[start]->fragment_string) + 1);
            strcpy(reps[count], frag_arr[start]->fragment_string);
            count++;
            free(cycle_nodes);
            continue;
        }
 
        // find the weakest edge on this cycle (cycle closing edge)
        // that edge will be opened, removed to form a chain
        int min_ov = INT_MAX;
        int min_idx = 0; // index into cycle_nodes[] whose outgoing edge is weakest
 
        for (int i = 0; i < cycle_len; i++) {
            int from = cycle_nodes[i];
            int to = cycle_nodes[(i + 1) % cycle_len];
            int ov = selected[from][to];
            if (ov < min_ov) {
                min_ov  = ov;
                min_idx = i;
            }
        }
 
        /*
         * chain starts at the node after the removed edge's head
         * i.e. if we remove edge (cycle_nodes[min_idx] -> cycle_nodes[min_idx+1])
         * the chain starts at cycle_nodes[min_idx+1]
         */
        int chain_start = (min_idx + 1) % cycle_len;
 
        // calculate total length of the representative string.
        size_t total_len = strlen(frag_arr[cycle_nodes[chain_start]]->fragment_string);
        for (int step = 1; step < cycle_len; step++) {
            int prev = cycle_nodes[(chain_start + step - 1) % cycle_len];
            int cur2 = cycle_nodes[(chain_start + step) % cycle_len];
            int ov = selected[prev][cur2];
            total_len += strlen(frag_arr[cur2]->fragment_string) - ov;
        }
 
        // build the representative string by walking the chain.
        char *rep = malloc(total_len + 1);
        assert(rep != NULL);
        strcpy(rep, frag_arr[cycle_nodes[chain_start]]->fragment_string);
 
        for (int step = 1; step < cycle_len; step++) {
            int prev = cycle_nodes[(chain_start + step - 1) % cycle_len];
            int cur2 = cycle_nodes[(chain_start + step) % cycle_len];
            int ov = selected[prev][cur2];
            strcat(rep, frag_arr[cur2]->fragment_string + ov);
        }
 
        reps[count] = rep;
        count++;
        free(cycle_nodes);
    }
 
    free(visited);
    *rep_count = count;
    return reps;
}

/* Section 4 - TGREEDY, greedy merge of representative strings
*
*  an overlap matrix on the representatives and run the standard
*  greedy with cycle avoidance to merge them into a single string
*/

// computes the overlap between suffix a and prefix of b
// returns the length of the longest such overlap

static int
compute_overlap(const char *a, const char *b)
{
    size_t la = strlen(a);
    size_t lb = strlen(b);
    size_t max_ov = la < lb ? la : lb;
 
    for (size_t k = max_ov; k > 0; k--) {
        if (strncmp(a + la - k, b, k) == 0)
            return (int)k;
    }
    return 0;
}
 
// DFS reachability check for cycle detection in the greedy path builder
static int
can_reach(int **graph, int n, int from, int to, int *vis)
{
    if (from == to) return 1;
    vis[from] = 1;
    for (int next = 0; next < n; next++) {
        if (graph[from][next] > 0 && !vis[next]) {
            if (can_reach(graph, n, next, to, vis))
                return 1;
        }
    }
    return 0;
}
 
// standard greedy with cycle avoidance on a set of strings
// returns a heap allocated single merged string
static char *
greedy_merge(char **strings, int count)
{
    if (count == 0)
        return strdup("");
    if (count == 1)
        return strdup(strings[0]);
 
    // build overlap matrix for the representatives
    int **ov = malloc(count * sizeof(int *));
    assert(ov != NULL);
    for (int i = 0; i < count; i++) {
        ov[i] = calloc(count, sizeof(int));
        assert(ov[i] != NULL);
    }
    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++)
            if (i != j)
                ov[i][j] = compute_overlap(strings[i], strings[j]);
 
    // collect and sort edges
    int edge_cap = count * count;
    Edge *edges = malloc(edge_cap * sizeof(Edge));
    int edge_count = 0;
    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++)
            if (i != j && ov[i][j] > 0) {
                edges[edge_count].from   = i;
                edges[edge_count].to     = j;
                edges[edge_count].weight = ov[i][j];
                edge_count++;
            }
 
    qsort(edges, edge_count, sizeof(Edge), edge_cmp_desc);
 
    // greedy selection with cycle avoidance
    int **sel = malloc(count * sizeof(int *));
    for (int i = 0; i < count; i++)
        sel[i] = calloc(count, sizeof(int));
 
    int *has_out = calloc(count, sizeof(int));
    int *has_in = calloc(count, sizeof(int));
    int *vis = malloc(count * sizeof(int));
    int added = 0;
 
    for (int e = 0; e < edge_count && added < count - 1; e++) {
        int u = edges[e].from;
        int v = edges[e].to;
        if (has_out[u] || has_in[v])
            continue;
 
        // cycle check - would adding u -> v close a cycle?
        for (int k = 0; k < count; k++) vis[k] = 0;
        if (can_reach(sel, count, v, u, vis))
            continue;
 
        sel[u][v] = edges[e].weight;
        has_out[u] = 1;
        has_in[v] = 1;
        added++;
    }
 
    // identify chain starts, nodes with no incoming selected edge
    // walk each chain and concatenate
    size_t total_len = 0;
    for (int i = 0; i < count; i++)
        total_len += strlen(strings[i]);
    // subtract overlaps
    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++)
            total_len -= sel[i][j];
 
    char *result = malloc(total_len + 1);
    assert(result != NULL);
    result[0] = '\0';
 
    int *chain_done = calloc(count, sizeof(int));
 
    for (int start = 0; start < count; start++) {
        if (has_in[start] || chain_done[start])
            continue;
 
        // walk the chain from this start node and append if no overlap between disconnected chains
        if (result[0] != '\0') {
        }
        strcat(result, strings[start]);
        chain_done[start] = 1;
 
        int cur = start;
        while (1) {
            int next = -1;
            for (int j = 0; j < count; j++) {
                if (sel[cur][j] > 0) {
                    next = j;
                    break;
                }
            }
            if (next == -1) break;
            strcat(result, strings[next] + sel[cur][next]);
            chain_done[next] = 1;
            cur = next;
        }
    }
 
    // handle any isolated nodes not yet included
    for (int i = 0; i < count; i++) {
        if (!chain_done[i])
            strcat(result, strings[i]);
    }
 
    // cleanup
    for (int i = 0; i < count; i++) {
        free(ov[i]);
        free(sel[i]);
    }
    free(ov);
    free(sel);
    free(edges);
    free(has_out);
    free(has_in);
    free(vis);
    free(chain_done);
 
    return result;
}

/* Section 5 - exact DP with bitmask for small n <= 18
*  
*  Held-Karp style DP over all subsets of fragments
*  dp[mask][i] = shortest string length covering the fragments in mask, ending with fragment i
*
*  Complexity: O(2^n * n^2) time, O(2^n * n) space
*  only called when n <= 18
*/

#define MAX_EXACT_N 18

// returns the optimal string via heap allocated bitmask DP
// frags is an array of n strings, ov_matrix is the n*n overlap matrix

static char *
dp_exact(char **frags, int n, int **ov_matrix)
{
    int full_mask = (1 << n) - 1;
    long total_states = (long)(1 << n) * n;
 
    // dp[mask][i] = minimum total length of superstring covering
    // the fragments in mask with fragment i placed last
    int **dp = malloc((1 << n) * sizeof(int *));
    int **parent_mask = malloc((1 << n) * sizeof(int *));
    int **parent_node = malloc((1 << n) * sizeof(int *));
    assert(dp != NULL && parent_mask != NULL && parent_node != NULL);
 
    for (int mask = 0; mask <= full_mask; mask++) {
        dp[mask] = malloc(n * sizeof(int));
        parent_mask[mask] = malloc(n * sizeof(int));
        parent_node[mask] = malloc(n * sizeof(int));
        assert(dp[mask] != NULL);
        for (int i = 0; i < n; i++) {
            dp[mask][i] = INT_MAX;
            parent_mask[mask][i] = -1;
            parent_node[mask][i] = -1;
        }
    }
 
    // base case - each single fragment
    for (int i = 0; i < n; i++) {
        dp[1 << i][i] = (int)strlen(frags[i]);
    }
 
    // fill DP
    for (int mask = 1; mask <= full_mask; mask++) {
        for (int last = 0; last < n; last++) {
            if (!(mask & (1 << last)))
                continue;
            if (dp[mask][last] == INT_MAX)
                continue;
 
            // try adding fragment next that is not yet in mask
            for (int next = 0; next < n; next++) {
                if (mask & (1 << next))
                    continue;
 
                int new_mask = mask | (1 << next);
                int new_len = dp[mask][last] + (int)strlen(frags[next]) - ov_matrix[last][next];
 
                if (new_len < dp[new_mask][next]) {
                    dp[new_mask][next] = new_len;
                    parent_mask[new_mask][next] = mask;
                    parent_node[new_mask][next] = last;
                }
            }
        }
    }
 
    // find the best ending node
    int best_len = INT_MAX;
    int best_last = 0;
    for (int i = 0; i < n; i++) {
        if (dp[full_mask][i] < best_len) {
            best_len = dp[full_mask][i];
            best_last = i;
        }
    }
 
    // backtrack to recover the permutation order
    int *order = malloc(n * sizeof(int));
    assert(order != NULL);
    int mask = full_mask;
    int node = best_last;
    for (int pos = n - 1; pos >= 0; pos--) {
        order[pos] = node;
        int pmask = parent_mask[mask][node];
        int pnode = parent_node[mask][node];
        mask = pmask;
        node = pnode;
    }
 
    // build the result string from the permutation
    char *result = malloc(best_len + 1);
    assert(result != NULL);
    strcpy(result, frags[order[0]]);
 
    for (int pos = 1; pos < n; pos++) {
        int prev = order[pos - 1];
        int cur = order[pos];
        int ov = ov_matrix[prev][cur];
        strcat(result, frags[cur] + ov);
    }
 
    // cleanup
    for (int mask2 = 0; mask2 <= full_mask; mask2++) {
        free(dp[mask2]);
        free(parent_mask[mask2]);
        free(parent_node[mask2]);
    }
    free(dp);
    free(parent_mask);
    free(parent_node);
    free(order);
 
    return result;
}

// Section 6 - main

int
main(int argc, char *argv[])
{
    if (argc != 2) {
        fprintf(stderr,
            "Usage: %s [ <input_file> | - ]\n"
            "  Reads fragments (one per line) and outputs the shortest superstring.\n",
            argv[0]);
        exit(1);
    }

    // phase 1 - read and preprocess common.c

    Fragment *frags = read_all_fragments(argv[1]);
    if (frags == NULL) {
        fprintf(stderr, "Error: no fragments read.\n");
        return 1;
    }
 
    frags = remove_substring_fragments(frags);
 
    int n;
    int **overlap = build_overlap_matrix(frags, &n);
    Fragment **frag_arr = fragments_to_array(frags, n);
 
    if (n == 0) {
        print_solution("");
        free(frag_arr);
        free_all_fragments(frags);
        return 0;
    }
 
    if (n == 1) {
        print_solution(frag_arr[0]->fragment_string);
        free(overlap[0]);
        free(overlap);
        free(frag_arr);
        free_all_fragments(frags);
        return 0;
    }

    // prepare a plain string array
        char **frag_strings = malloc(n * sizeof(char *));
    for (int i = 0; i < n; i++)
        frag_strings[i] = frag_arr[i]->fragment_string;
 
    char *best_solution = NULL;

    // phase 2-4 - MGREEDY + open cycles + TGREEDY
    {
        int **cc = mgreedy_cycle_cover(overlap, n);
 
        int rep_count = 0;
        char **reps = open_cycles_to_representatives(cc, n, frag_arr, &rep_count);
 
        char *tgreedy_result = greedy_merge(reps, rep_count);
 
        if (solution_is_valid(tgreedy_result, frags)) {
            best_solution = tgreedy_result;
        } else {
            // safety fallback - if TGREEDY somehow produces an invalid result, fall back to naive concatenation.
            fprintf(stderr, "TGREEDY result invalid, falling back.\n");
            free(tgreedy_result);
        }
 
        for (int i = 0; i < rep_count; i++)
            free(reps[i]);
        free(reps);
        for (int i = 0; i < n; i++)
            free(cc[i]);
        free(cc);
    }

    // phase 5 - DP for small inputs (Held-Karp algorithm)
    if (n <= MAX_EXACT_N) {
        char *exact_result = dp_exact(frag_strings, n, overlap);
 
        if (solution_is_valid(exact_result, frags)) {
            if (best_solution == NULL ||
                strlen(exact_result) < strlen(best_solution)) {
                free(best_solution);
                best_solution = exact_result;
            } else {
                free(exact_result);
            }
        } else {
            fprintf(stderr, "DP result invalid. Cleaning up\n");
            free(exact_result);
        }
    }

    // fallback - if everything failed, just do naive concatenation
    if (best_solution == NULL) {
        size_t total = 0;
        for (int i = 0; i < n; i++)
            total += strlen(frag_strings[i]);
        best_solution = malloc(total + 1);
        best_solution[0] = '\0';
        for (int i = 0; i < n; i++)
            strcat(best_solution, frag_strings[i]);
    }

    // output
    print_solution(best_solution);

    // cleanup
    free(best_solution);
    free(frag_strings);
    for (int i = 0; i < n; i++)
        free(overlap[i]);
    free(overlap);
    free(frag_arr);
    free_all_fragments(frags);
 
    return 0;
}
