#ifndef SEMANT_H_
#define SEMANT_H_

#include <assert.h>
#include <iostream>  
#include <map>
#include <set>
#include <vector>
#include <utility>
#include "cool-tree.h"
#include "stringtab.h"
#include "symtab.h"
#include "list.h"

#define TRUE 1
#define FALSE 0


// This is a structure that may be used to contain the semantic
// information such as the inheritance graph.  You may use it or not as
// you like: it is only here to provide a container for the supplied
// methods.

struct ClassTreeNode {
	bool inited;
	Symbol parent;
	Class_ class_;
	std::vector<Symbol> children;

	ClassTreeNode() 
	{
		inited = false;
	}
};


#endif

