//
// The following include files must come first.

#ifndef COOL_TREE_HANDCODE_H
#define COOL_TREE_HANDCODE_H

#include <iostream>
#include <map>
#include <vector>
#include <utility>
#include "tree.h"
#include "cool.h"
#include "stringtab.h"
#define yylineno curr_lineno;
extern int yylineno;

inline Boolean copy_Boolean(Boolean b) {return b; }
inline void assert_Boolean(Boolean) {}
inline void dump_Boolean(ostream& stream, int padding, Boolean b)
	{ stream << pad(padding) << (int) b << "\n"; }

void dump_Symbol(ostream& stream, int padding, Symbol b);
void assert_Symbol(Symbol b);
Symbol copy_Symbol(Symbol b);

class Program_class;
typedef Program_class *Program;
class Class__class;
typedef Class__class *Class_;
class Feature_class;
typedef Feature_class *Feature;
class Formal_class;
typedef Formal_class *Formal;
class Expression_class;
typedef Expression_class *Expression;
class Case_class;
typedef Case_class *Case;

typedef list_node<Class_> Classes_class;
typedef Classes_class *Classes;
typedef list_node<Feature> Features_class;
typedef Features_class *Features;
typedef list_node<Formal> Formals_class;
typedef Formals_class *Formals;
typedef list_node<Expression> Expressions_class;
typedef Expressions_class *Expressions;
typedef list_node<Case> Cases_class;
typedef Cases_class *Cases;

#define Program_EXTRAS                          \
virtual void cgen(ostream&) = 0;		\
virtual void dump_with_types(ostream&, int) = 0; 



#define program_EXTRAS                          \
void cgen(ostream&);     			\
void dump_with_types(ostream&, int);            

#define Class__EXTRAS                   \
virtual Symbol get_name() = 0;  	\
virtual Symbol get_parent() = 0;    	\
virtual Symbol get_filename() = 0;      \
virtual void dump_with_types(ostream&,int) = 0; 


#define class__EXTRAS                                  \
Symbol get_name()   { return name; }		       \
Symbol get_parent() { return parent; }     	       \
Symbol get_filename() { return filename; }             \
void dump_with_types(ostream&,int);                    


#define Feature_EXTRAS                                        \
virtual void dump_with_types(ostream&,int) = 0; \
virtual int get_type() = 0; \
virtual void set_offset(int, Symbol &) = 0; \
virtual void initialization(ostream &) { } \
virtual void emit_method(ostream &) {} \
virtual Symbol get_name() = 0; \
int get_offset() { return offset; } \
int offset; \
static const int ATTR = 0; \
static const int METHOD = 1;


#define Feature_SHARED_EXTRAS                                       \
void dump_with_types(ostream&,int);   \
int get_type(); \
Symbol get_name() { return name; }\
void set_offset(int, Symbol &);

#define attr_EXTRAS \
void initialization(ostream &s);

#define method_EXTRAS \
void emit_method(ostream &);



#define Formal_EXTRAS                              \
virtual void dump_with_types(ostream&,int) = 0;\
virtual Symbol get_name() = 0;


#define formal_EXTRAS                           \
void dump_with_types(ostream&,int);\
Symbol get_name() { return name; }


#define Case_EXTRAS                             \
virtual void dump_with_types(ostream& ,int) = 0;


#define branch_EXTRAS                                   \
void dump_with_types(ostream& ,int);


#define Expression_EXTRAS                    \
Symbol type;                                 \
Symbol get_type() { return type; }           \
Expression set_type(Symbol s) { type = s; return this; } \
virtual void code(ostream&) = 0; \
virtual bool is_no_expr() { return false; }; \
virtual void dump_with_types(ostream&,int) = 0;  \
void dump_type(ostream&, int);               \
Expression_class() { type = (Symbol) NULL; }

#define Expression_SHARED_EXTRAS           \
void code(ostream&); 			   \
void dump_with_types(ostream&,int); 

#define no_expr_EXTRAS \
bool is_no_expr() { return true; }



#endif
