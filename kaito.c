/* compile with:
    gcc -O2 -o graph_reconstruction graph_text_reconstruction.c
*/

/* Overall Complexity: O(n^4 + n^2 * m^2)
    n = number of input fragments
    m = maximum length of any single fragment
    L = total length of all fragments (<= n * m)
   Dominant phases:
    select_overlap_edges        : O(n^4)
    build_overlap_matrix        : O(n^2 * m^2)
    remove_substring_fragments  : O(n^2 * m^2)
*/

#include <stdio.h>     // getline(), fputs(), fputc(), fprintf(), printf(), fopen(), fclose(), FILE, stdin, stdout, stderr
#include <stdlib.h>    // malloc(), calloc(), free(), qsort(), exit()
#include <string.h>    // strcmp(), strncmp(), strstr(), strlen(), strcpy(), strcat()
#include <assert.h>    // assert()
#include <stdbool.h>   // bool

struct fragment_s {            // datatype "struct fragment_s" for the elements of a singly linked list.
  struct fragment_s * next_fragment;
  char * fragment_string;
};

typedef struct {
  int i;
  int j;
  int weight;
} Edge;

static int
compare_desc(const void *a, const void *b)
{
  const Edge *ea = (const Edge *)a;
  const Edge *eb = (const Edge *)b;
  return eb->weight - ea->weight;
}

static int
reachable(int **graph, int n, int from, int to, int *visited)
{
  if(from == to) return 1;
  visited[from] = 1;
  for(int next = 0; next < n; next++) {
    if(graph[from][next] > 0 && !visited[next]) {
      if (reachable(graph, n, next, to, visited)) return 1;
    }
  }
  return 0;
}

struct fragment_s *
read_all_fragments( char const * file_name )
{
  FILE * input = stdin;                                                   // read from standard input if file_name is "-",
  if( strcmp( file_name, "-" ) != 0 )
    input = fopen( file_name, "r" );                                       // otherwise open the named file for "r"eading.
  if( input == NULL ) {
    fprintf( stderr, "Error: file could not be opened for reading: %s\n", file_name );
    return NULL;
  }
  struct fragment_s * top_fragment = NULL;    // prepare a singly linked list to hold the text fragments, initially empty.
  while( 1 ) {
    char * line_buffer = NULL;           // memory buffer for the next line read from input, to be allocated by getline().
    size_t buffer_size = 0;            // this variable will be set to the size of the line buffer allocated by getline().
    ssize_t read_bytes = getline( & line_buffer, & buffer_size, input );     // allocate buffer and store next input line.
    if( read_bytes <= 0 ) {
      if( line_buffer != NULL ) free( line_buffer );  // 'man getline' says the caller must free the buffer even on error.
      break;           // break from the input read loop as soon as nothing else can be read, such as at the end of input.
    }
    assert( line_buffer != NULL );    // ensure getline() did allocate a fresh buffer and save its address in line_buffer.
    assert( line_buffer[read_bytes] == '\0' );  // ensure the getline()-allocated string ends with extra '\0' as promised.
    int i = 0;
    for( i = 0; i < read_bytes; i ++ )  // some sanitising required: since getline() saves any read newline in the buffer,
      if( line_buffer[i] == '\n' || line_buffer[i] == '\r' )       // we replace it (on Linux/MacOS/Windows) with an '\0'.
        line_buffer[i] = '\0';
    struct fragment_s * new_fragment = malloc( sizeof( struct fragment_s ));
    new_fragment->next_fragment = top_fragment;
    new_fragment->fragment_string = line_buffer;                                                   // passing the baton...
    line_buffer = NULL;                                                                           // don't keep the baton!
    top_fragment = new_fragment;
#ifdef DEBUG
printf( "%ld: (%ld) %s\n", read_bytes, buffer_size, top_fragment->fragment_string );
#endif
  }
  fclose( input );
  input = NULL;
  return top_fragment;
}

void
free_all_fragments( struct fragment_s * top_fragment )
{
  while( top_fragment != NULL ) {
    struct fragment_s * this_fragment = top_fragment;
    top_fragment = this_fragment->next_fragment;
    free( this_fragment->fragment_string );
    free( this_fragment);
  }
}

struct fragment_s *
remove_substring_fragments(struct fragment_s * top_fragment)
{
  struct fragment_s * prev_fragment = NULL;
  struct fragment_s * test_fragment = top_fragment;
  while(test_fragment != NULL) {
    int removed = 0;
    struct fragment_s * search_fragment = top_fragment;
    while(search_fragment != NULL) {
      if(search_fragment != test_fragment && strstr(search_fragment->fragment_string, test_fragment->fragment_string) !=NULL) {
        struct fragment_s * to_free = test_fragment;
        if(prev_fragment == NULL) {
          top_fragment = test_fragment->next_fragment;
          test_fragment = top_fragment;
        }
        else {
          prev_fragment->next_fragment = test_fragment->next_fragment;
          test_fragment = test_fragment->next_fragment;
        }
        free(to_free->fragment_string);
        free(to_free);
        removed = 1;
        break;
      }
      search_fragment = search_fragment->next_fragment;
    }
    if(!removed) {
      prev_fragment = test_fragment;
      test_fragment = test_fragment->next_fragment;
    }
  }
  return top_fragment;
}

int **
build_overlap_matrix(struct fragment_s * top_fragment, int * fragment_count)
{
  int n = 0;
  for(struct fragment_s * p = top_fragment; p != NULL; p = p->next_fragment) n++;

  int ** graph = malloc(n * sizeof(int *));
  for(int i = 0; i < n; i++) {
    graph[i] = calloc(n, sizeof(int));
  }

  int i = 0;
  for(struct fragment_s * test_fragment = top_fragment; test_fragment != NULL; test_fragment = test_fragment->next_fragment) {
    char * test_string = test_fragment->fragment_string;
    size_t test_length = strlen(test_string);
    int j = 0;
    for(struct fragment_s * search_fragment = top_fragment; search_fragment != NULL; search_fragment = search_fragment->next_fragment) {
      if(test_fragment != search_fragment) {
        char * search_string = search_fragment->fragment_string;
        size_t search_length = strlen(search_string);
        size_t min_length = test_length < search_length ? test_length : search_length;
        for(size_t k = min_length; k > 0; k--) {
          if(strncmp(test_string + test_length - k, search_string, k) == 0) {
            graph[i][j] = (int)k;
            break;
          }
        }
      }
      j++;
    }
    i++;
  }
  *fragment_count = n;
  return graph;
}

int **
select_overlap_edges(int **overlap_matrix, int fragment_count)
{
  int edge_count = 0;
  Edge *edges = malloc(fragment_count * fragment_count * sizeof(Edge));
  for(int i = 0; i < fragment_count; i++) {
    for(int j = 0; j < fragment_count; j++) {
      if(i != j && overlap_matrix[i][j] > 0) {
        edges[edge_count].i = i;
        edges[edge_count].j = j;
        edges[edge_count].weight = overlap_matrix[i][j];
        edge_count++;
      }
    }
  }

  qsort(edges, edge_count, sizeof(Edge), compare_desc);

  int **selected_graph = malloc(fragment_count * sizeof(int *));
  for(int i = 0; i < fragment_count; i++) {
    selected_graph[i] = calloc(fragment_count, sizeof(int));
  }
  
  int edge_added = 0;
  int edge_index = 0;
  int *visited = malloc(fragment_count * sizeof(int));
  while(edge_added < fragment_count - 1 && edge_index < edge_count) {
    Edge test_edge = edges[edge_index];
    int addable = 1;

    for(int k = 0; k < fragment_count; k++) {
      if (selected_graph[test_edge.i][k] > 0) { addable = 0; break; }
      if (selected_graph[k][test_edge.j] > 0) { addable = 0; break; }
    }

    if(addable) {
      for(int k = 0; k < fragment_count; k++) visited[k] = 0;
      if(reachable(selected_graph, fragment_count, test_edge.j, test_edge.i, visited)) {
        addable = 0;
      }
    }

    if(addable) {
      selected_graph[test_edge.i][test_edge.j] = test_edge.weight;
      edge_added++;
    }
    edge_index++;
  }

  free(visited);
  free(edges);
  return selected_graph;
}

char **
build_sequence_list(struct fragment_s * top_fragment, int ** overlap_matrix_selected, int fragment_count, size_t * sequence_count)
{
  char ** sequence_list = malloc(fragment_count * sizeof(char *));
  size_t count = 0;

  struct fragment_s ** fragment_array = malloc(fragment_count * sizeof(struct fragment_s *));
  struct fragment_s * this_fragment = top_fragment;
  for(int i = 0; i < fragment_count; i++) {
    fragment_array[i] = this_fragment;
    this_fragment = this_fragment->next_fragment;
  }

  for(int j = 0; j < fragment_count; j++) {
    int is_start = 1;
    for(int k = 0; k < fragment_count; k++) {
      if(overlap_matrix_selected[k][j] > 0) {
        is_start = 0;
        break;
      }
    }
    if(!is_start) continue;

    size_t total_len = strlen(fragment_array[j]->fragment_string);
    {
      int from = j;
      while(1) {
        int next = -1;
        for(int l = 0; l < fragment_count; l++) {
          if(overlap_matrix_selected[from][l] > 0) {
            next = l;
            break;
          }
        }
        if(next == -1) break;
        total_len += strlen(fragment_array[next]->fragment_string)
                   - overlap_matrix_selected[from][next];
        from = next;
      }
    }

    char * sequence = malloc(total_len + 1);
    strcpy(sequence, fragment_array[j]->fragment_string);

    int from = j;
    while(1) {
      int next = -1;
      for(int l = 0; l < fragment_count; l++) {
        if(overlap_matrix_selected[from][l] > 0) {
          next = l;
          break;
        }
      }
      if(next == -1) break;

      int overlap = overlap_matrix_selected[from][next];
      strcat(sequence, fragment_array[next]->fragment_string + overlap);
      from = next;
    }

    sequence_list[count] = malloc(strlen(sequence) + 1);
    strcpy(sequence_list[count], sequence);
    count++;

    free(sequence);
  }

  free(fragment_array);

  *sequence_count = count;
  return sequence_list;
}

int
main(int argc, char * argv[])
{
  // Validate command-line arguments.
  // Complexity: O(1)
  if(argc != 2) {
    fprintf( stderr, "Usage: %s [ <input_file> | - ]\n with input_file (or standard input if -) containing all string fragments on successive lines (no blanks)\n", argv[0] );
    exit(1);
  }

  // Build a linked list of fragments read from input.
  // Complexity: O(L)
  struct fragment_s * top_fragment = read_all_fragments(argv[1]);

  // Remove any fragment that is already a substring of another fragment.
  // Complexity: O(n^2 * m^2)
  top_fragment = remove_substring_fragments(top_fragment);

  // Represent fragments as vertices in a weighted directed graph, with edge weights equal to the number of overlapping characters between fragments.
  // Complexity: O(n^2 * m^2)
  int fragment_count;
  int ** overlap_matrix = build_overlap_matrix(top_fragment, &fragment_count);

  // Select edges by descending weight, keeping each vertex's in/out-degree ≤ 1 and avoiding cycles; discard the rest.
  // Complexity: O(n^4)
  int ** overlap_matrix_selected = select_overlap_edges(overlap_matrix, fragment_count);

  // For each component of the weighted directed graph, traverse in topological order to build an ordered, overlap-free sequence of fragments.
  // Complexity: O(n^2 * m)
  size_t sequence_count;
  char ** sequence_list = build_sequence_list(top_fragment, overlap_matrix_selected, fragment_count, &sequence_count);

  // Merge each sequence into a single string.
  // Complexity: O(L^2)
  int merged_len = 0;
  for(int i = 0; i < sequence_count; i++) {
    merged_len += strlen(sequence_list[i]);
  }

  char * merged = malloc(merged_len + 1);
  merged[0] = '\0';
  for(int i = 0; i < sequence_count; i++) {
    strcat(merged, sequence_list[i]);
  }
  fputs(merged, stdout);
  fputc('\n', stdout);

  // Cleanup
  // Complexity: O(n^2)
  for(int i = 0; i < fragment_count; i++) {
    free(overlap_matrix[i]);
    free(overlap_matrix_selected[i]);
  }
  free(overlap_matrix);
  free(overlap_matrix_selected);
  for(int i = 0; i < sequence_count; i++) {
    free(sequence_list[i]);
  }
  free(sequence_list);
  free_all_fragments(top_fragment);
  free(merged);
  return 0;
}