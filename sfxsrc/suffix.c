// Author: Patrick Brodie


#include "suffix.h"

// ============================================================================
// suffix.c contains the implementation of McCreight's linear-time suffix tree 
// building algorithm.  
// ============================================================================


void allocate_node (struct node **node, int sufnum, int starti, 
					int endi, struct node *parent)
// Allocate one node which marks the edge corresponding to the slice 
// input_string[starti: endi].  
{
	*node = (struct node*) malloc (sizeof (struct node));
	(*node) -> id = idCnt++;
	(*node) -> sfxnum = sufnum;
	(*node) -> strdepth = parent -> strdepth + endi - starti;
	(*node) -> starti = starti;
	(*node) -> endi = endi;
	(*node) -> array_start = -1;
	(*node) -> array_end = -1;
	(*node) -> sfxlink = NULL;
	(*node) -> leftchild = NULL;
	(*node) -> rightsib = NULL;
	(*node) -> parent = parent;

	// Check exact matching sequence length
	if (!deepest) {
		deepest = root;
	} else if (sufnum == -1 && (*node) -> strdepth > deepest -> strdepth) {
		deepest = (*node);
	}
}


void init_root (int length)
// Initialize the root node of the suffix tree.
{
	// Give root node unique (no edge) attributes
	root = (struct node*) malloc (sizeof(struct node));
	root -> id = 0;
	root -> sfxnum = -1;
	root -> strdepth = 0;
	root -> starti = -1;
	root -> endi = -1;

	// suffix link of root is itself.
	root -> sfxlink = root;

	// Initialize the first child of root with suffix corresponding to the
	// entire input string.
	allocate_node (&(root -> leftchild), 0, 0, length, root);

	// Root has no sibs or parent.
	root -> rightsib = NULL;
	root -> parent = NULL;
}


void prepare_str (char *s) 
// Prepare the string to be inserted into the tree.
// Store its length and point input_string at it.
{
	input_string = s;
	if (s[strlen(s) - 1] == '$') {	// User has already prepared string.
		slen = strlen (s) - 1;
	} else {						// String has not yet been prepared.
		slen = strlen (s);
		input_string [slen] = '$';
		input_string [slen + 1] = 0;
	}

}


void print_string_slice (char *s, int start, int end)
// Print the slice of string s[start:end]
{
	while (start < end) {
		printf ("%c", s[start++]);
	}
	printf ("\n");
}


void print_sequence (char *name, char *sequence)
// Output the given sequence with 50 chars per line.
{
	int cnt = 0;
	printf ("========= Printing Sequence: "); printf ("%s\n\n", name);
	while (*sequence) {
		if (cnt == 50) {
			printf ("\n");
			cnt = 0;
		}
		printf ("%c", *sequence++);
		++cnt;
	}
	printf ("\n");
}


void print_children (struct node *node)
// Print the children of the given node from left to right.
{
	struct node *child;
	child = node -> leftchild;
	printf ("Children of node %d\n", node -> id);
	while (child) {
		printf ("ID: %d | Depth: %d | Interval: [%d, %d)\n", 
			child -> id, child -> strdepth, child -> starti, child -> endi);
		child = child -> rightsib;
	}
}


void print_BWT (struct node *tree)
// Print the BWT index for the input string using the constructed tree.
{
	struct node *curr;
	if (tree) {
		curr = tree -> leftchild;
		while (curr) {
			print_BWT (curr);
			curr = curr -> rightsib;
		}
		if (tree -> sfxnum != -1) {
			// BWT index is '$' if sfxnum == 0, else input_string[sfxnum - 1]
			if (tree -> sfxnum == 0) {
				printf ("$\n");
			} else {
				printf ("%c\n", input_string[tree -> sfxnum - 1]);
			}
		}
	}
}


void print_DFS (struct node *tree)
// Print string-depth info for depth-first traversal of the given suffix tree.
{
	struct node *curr;
	if (tree) {
		curr = tree -> leftchild;
		while (curr) {
			print_DFS (curr);
			curr = curr -> rightsib;
		}
		printf ("%4d",tree -> strdepth);
	}
}


void do_DFS (struct node *tree) 
// Perform a depth-first traversal of the given suffix tree
// and store the string depths of the nodes in the order visited.
{
	struct node *curr;
	if (tree) {
		curr = tree -> leftchild;
		while (curr) {
			do_DFS (curr);
			curr = curr -> rightsib;
		}
		if (tree -> sfxnum == -1) {
			++numints;
		} else {
			++numleaves;
		}
	}
}


struct node *find_node (struct node *tree, int sufnum) 
// Locate the given node within the tree and return a pointer to it.
{
	struct node *curr, *ret;
	if (!tree) {
		return NULL;
	} else if (tree -> sfxnum == sufnum) {
		return tree;
	} else {
		curr = tree -> leftchild;
		while (curr != NULL) {
			if (ret = find_node (curr, sufnum))
				return ret;
			curr = curr -> rightsib;
		}
		return NULL;
	}
}


void push_node (struct node **list, struct node *node)
// Push the given node onto the front of the given list.
{
	node -> rightsib = *list;
	*list = node;
}


void sorted_insert (struct node **list, struct node *node)
// Insert the given node into the given list in alphabetical order according
// to input_string[node -> starti].
{
	struct node **curr = list;
	if (!(input_string[node -> starti] == '$')) {	// put $ at the start of list
		while (*curr && input_string[(*curr) -> starti] < input_string[node -> starti]) {
			curr = &((*curr) -> rightsib);
		}
	}
	push_node (curr, node);
}


void push_branch (struct node **tree, struct node *new_child)
// Push the given new child to create a new branch off the given tree.
{
	sorted_insert (&((*tree) -> leftchild), new_child);
}


void push_int_node (struct node **tree, struct node *node)
// Push a node into the given tree.
// This is used to "break" edges to create new internal nodes.
{
	node -> leftchild = (*tree) -> leftchild;
	(*tree) -> leftchild = node;
	(*tree) -> endi = node -> starti;
}


struct node *remove_child (struct node **list, struct node *child)
// Remove a child from the given list of children and return a pointer to it.
{
	struct node *tmp, **curr;
	curr = list;
	while (*curr) {
		if ((*curr) -> id == child -> id) {
			break;
		}
		curr = &((*curr) -> rightsib);
	}
	tmp = *curr;
	*curr = (*curr) -> rightsib;
	if (tmp) {
		tmp -> rightsib = NULL;
		tmp -> parent = NULL;
	}
	return tmp;
}


int identify_instype (struct node *u)
// Identify the insertion type for the parent u of the last inserted leaf.
// 		IA: suffix link of u is known and u is not root.
// 		IB: suffix link of u is known and u is root.
//		IIA: suffix link of u is unknown and u's parent is not root.
// 		IIB: suffix link of u is unknown and u's parent is root.
{
	if (u -> sfxlink) {
		return (u == root)? IB : IA;
	} else {
		return (u -> parent == root)? IIB : IIA;
	}
}


struct node *get_branch_by_match (char c, struct node *parent)
// Search the children of the given parent node for the child whose label
// begins with the given character c.  Return it when found, NULL if not found.
{
	struct node *branch = parent -> leftchild;
	while (branch) {
		if (c == input_string[branch -> starti]) {
			break;
		}
		branch = branch -> rightsib;
	}
	return branch;
}


struct node *get_branch (int matchindex, struct node *parent)
// Search the children of the given parent node for the child whose label
// begins with input_string[matchindex].  Return it when found, NULL if not found.
{
	struct node *branch = parent -> leftchild;
	while (branch) {
		if (input_string[matchindex] == input_string[branch -> starti]) 
			break;
		branch = branch -> rightsib;
	}
	return branch;
}


struct node *break_edge (int breakindex, struct node *breaknode)
// Break an edge by inserting a new node labelled with the first portion of the broken edge,
// Push the rest of the edge label as a child branch of the new node.
// Return a reference to the new internal node.
{
	struct node *newint;
	allocate_node (&newint, -1, breaknode -> starti, breakindex, breaknode -> parent);
	breaknode = remove_child (&(breaknode -> parent -> leftchild), breaknode);
	breaknode -> parent = newint;
	breaknode -> starti = newint -> endi;
	push_branch (&newint, breaknode);
	push_branch (&(newint -> parent), newint);
	return newint;
}


struct node *find_path (int index, int sufdepth, struct node *v)
// Find path to the insertion point for the suffix beginning at input_string[index].
// Allocate and insert the node.
{
	int i, j, e;
	struct node *branch, *parent, *newint, *leaf, *tmp;
	
	// find child starting with input_string[sufdepth]
	i = sufdepth;
	parent = v;
	branch = get_branch (i, v);
	if (branch) {
		j = branch -> starti;
		e = branch -> endi - branch -> starti;
		// Look for the first mismatch in the current suffix and the path below v
		while (input_string[i] == input_string[j]) {
			if (!(--e)) {
				parent = branch;
				if ((branch = get_branch (++i, branch)) == NULL) break;
				j = branch -> starti;
				e = branch -> endi - branch -> starti;
			} else {
				++i; ++j;
			}
		}
	}
	
	// Allocate and insert a leaf at the branch point.
	if (branch) {	// Fork branch by making new internal node
		if (branch -> endi - branch -> starti > 1) {
			newint = break_edge (j, branch);
			allocate_node (&leaf, index, i, slen, newint);
			push_branch (&newint, leaf);
		} else {
			allocate_node (&leaf, index, i, slen, branch);
			push_branch (&branch, leaf);
		}
	} else {
		allocate_node (&leaf, index, i, slen, parent);
		push_branch (&parent, leaf);
	}
	return leaf;
}


struct node *consume_beta (int index, int firsti, int betalen, 
							struct node *startnode, struct node *u)
// Consume Beta (slice of input string);
// Establish suffix link of u to point to the spot at which Beta was consumed;
// Find leaf insertion location and return a pointer to it.
{
	int e, r, depth;
	struct node *branch, *leaf;

	r = 0;
	branch = startnode;

	while (r < betalen) {
		e = (branch -> endi) - (branch -> starti);
		if (r + e > betalen) {	// Beta is consumed mid-edge.
			// Break edge.  Assign the suffix link to the breakpoint and append leaf.
			u -> sfxlink = break_edge (branch -> starti + (betalen - r), branch); 
			allocate_node (&leaf, index, index + u -> sfxlink -> strdepth, slen, u -> sfxlink);
			push_branch (&(u -> sfxlink), leaf);
			break;
		} else if (r + e == betalen) {	// Beta is consumed at an existing node
			// Establish link and place the leaf by matching characters.
			u -> sfxlink = branch;
			depth = index - 1 + u -> strdepth;
			leaf = find_path (index, depth, branch);
			break;
		} else { 		// Hop to next node by adding its edge length to r.
			r += e;		
			branch = get_branch (firsti + r, branch);
		}
	}
	return leaf;
}

struct node *find_link (int index, struct node *parent, 
						struct node *child, int firsti, int lasti)
// Find and establish the suffix link for the given child node, 
// working downward from its parent's suffix link.
{
	int e, betalen, r, depth;
	struct node *vp, *branch, *leaf;

	// Capture beta for later consumption
	betalen = lasti - firsti;

	// traverse parent's suffix link.
	vp = parent -> sfxlink;
	branch = get_branch (firsti, vp);  

	if (branch && betalen) {	// Beta not empty
		// Consume Beta and Establish child -> sfxlink
		leaf = consume_beta (index, firsti, betalen, branch, child);
	} else {	
		// Beta was empty or no branch exists for the letter at input[index]
		child -> sfxlink = vp;
		leaf = find_path (index, index, vp);
	}
	return leaf;
}



// ---------------- Case Handlers ------------------------------

struct node *handle_IA (int index, struct node **tree, struct node *u) 
// Handle the case in which the suffix link of node u is KNOWN, and u is not root.
{
	int imin1, k, depth;
	struct node *v;
	k = u -> strdepth;
	imin1 = index - 1;
	v = u -> sfxlink;
	if (v == root) {
		depth = index;
	} else {
		depth = k + imin1;
	}
	// Find the path to the insertion point, insert, and return pointer to it.
	return find_path (index, depth, v);
}

struct node *handle_IB (int index, struct node **tree, struct node *u)
// Handle the case in which the suffix link of node u is KNOWN, and u is root.
{
	// Find the path to the insertion point, insert, and return pointer to it.
	return find_path (index, index, u);
}

struct node *handle_IIA (int index, struct node **tree, struct node *u)
// Handle the case in which the suffix link of node u is UNKNOWN, and u's parent is not root.
{
	struct node *up = u -> parent;
	// B = input_string [u->starti : u->endi]
	return find_link (index, up, u, u -> starti, u -> endi);

}

struct node *handle_IIB (int index, struct node **tree, struct node *u)
// Handle the case in which the suffix link of node u is UNKNOWN, and u's parent is root.
{
	struct node *up = u -> parent;
	return find_link (index, up, u, u -> starti + 1, u -> endi);
}

// =============================================================


struct node* insert_suffix (int index, struct node *tree, struct node **lastleaf)
// insert a new node corresponding to the suffix which begins at the given index.
{	
	struct node *u;
	int instype;

	u = (*lastleaf) -> parent;
	instype = identify_instype (u);
	switch (instype) {
		// Suffix Link of u is KNOWN and u is NOT root.
		case IA: 	*lastleaf = handle_IA (index, &tree, u); break;
		// Suffix Link of u is KNOWN and u is root.
		case IB: 	*lastleaf = handle_IB (index, &tree, u); break;
		// Suffix Link of u is UNKNOWN and u's parent is NOT root.
		case IIA:	*lastleaf = handle_IIA (index, &tree, u); break;
		// Suffix Link of u is UNKNOWN and u's parent is root.
		case IIB:	*lastleaf = handle_IIB (index, &tree, u); break;
		default: printf ("Something is rotten in the state of Denmark.\n"); break;
	}

	return tree;
}


struct node *build_tree (char *s, char *alphabet)
// Build a suffix tree for the given string over the given alphabet
// using McCreight's Suffix Link algorithm.
{
	int index = 1;
	struct node *lastleaf;
	idCnt = 1;
	
	prepare_str (s);
	if (s) {
		init_root (slen);
		if (slen > 0) {
			lastleaf = root -> leftchild;
			// Iterate over the input string s inserting each suffix into
			// the tree pointed to by root.
			while (s[index]) {
				root = insert_suffix (index++, root, &lastleaf);
			}
		}
	}
	return root;
}


void free_tree (struct node *tree) 
// Deallocate the memory allocated to the given tree structure
{
	struct node *tmp, *curr;
	if (tree) {
		curr = tree -> leftchild;
		while (curr) {
			tmp = curr;
			curr = curr -> rightsib;
			free_tree (tmp);
		}
		free (tree);
	}
}








