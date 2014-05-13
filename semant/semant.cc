

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <map>
#include <set>
#include <utility>
#include "semant.h"
#include "utilities.h"


extern int semant_debug;
extern char *curr_filename;

//////////////////////////////////////////////////////////////////////
//
// Symbols
//
// For convenience, a large number of symbols are predefined here.
// These symbols include the primitive type and method names, as well
// as fixed names used by the runtime system.
//
//////////////////////////////////////////////////////////////////////
static Symbol 
    arg,
    arg2,
    Bool,
    concat,
    cool_abort,
    copy,
    Int,
    in_int,
    in_string,
    IO,
    length,
    Main,
    main_meth,
    No_class,
    No_type,
    Object,
    out_int,
    out_string,
    prim_slot,
    self,
    SELF_TYPE,
    Str,
    str_field,
    substr,
    type_name,
    val;
//
// Initializing the predefined symbols.
//
static void initialize_constants(void)
{
    arg         = idtable.add_string("arg");
    arg2        = idtable.add_string("arg2");
    Bool        = idtable.add_string("Bool");
    concat      = idtable.add_string("concat");
    cool_abort  = idtable.add_string("abort");
    copy        = idtable.add_string("copy");
    Int         = idtable.add_string("Int");
    in_int      = idtable.add_string("in_int");
    in_string   = idtable.add_string("in_string");
    IO          = idtable.add_string("IO");
    length      = idtable.add_string("length");
    Main        = idtable.add_string("Main");
    main_meth   = idtable.add_string("main");
    //   _no_class is a symbol that can't be the name of any 
    //   user-defined class.
    No_class    = idtable.add_string("_no_class");
    No_type     = idtable.add_string("_no_type");
    Object      = idtable.add_string("Object");
    out_int     = idtable.add_string("out_int");
    out_string  = idtable.add_string("out_string");
    prim_slot   = idtable.add_string("_prim_slot");
    self        = idtable.add_string("self");
    SELF_TYPE   = idtable.add_string("SELF_TYPE");
    Str         = idtable.add_string("String");
    str_field   = idtable.add_string("_str_field");
    substr      = idtable.add_string("substr");
    type_name   = idtable.add_string("type_name");
    val         = idtable.add_string("_val");
}

static int semant_errors = 0;
static ostream& error_stream = cerr;
static Classes classes_;

static const int TYPE_CLASS  = 1;
static const int TYPE_SELF   = 2;
static const int TYPE_INSIDE = 3;

void install_basic_classes();
bool HasCycle(Symbol, std::map<Symbol, int> &);
void TraversalClass(Symbol);
void CreateMethodTable(Symbol, const std::vector<Feature> &);
int CheckTypeClass(Symbol);
bool TypeCheck(Symbol, Symbol);
Symbol Lub(Symbol, Symbol);
Symbol RegularType(Symbol, bool &);

Symbol SetExprType(Expression);
Symbol CommonDispatchExpr(Symbol, Symbol, Symbol, Expressions);
Symbol UnaryExpr(Expression, Symbol);
Symbol BinaryExpr(Expression, Expression, Symbol);
Symbol Fmt_SELF_TYPE(Symbol);

// own code
typedef std::map<Symbol, ClassTreeNode> ClassTree;

typedef std::pair<Symbol, Symbol> MethodTableKey;
typedef std::pair<Formals, Symbol> MethodTableValue;
typedef SymbolTable<MethodTableKey, MethodTableValue> MethodTable;

typedef SymbolTable<Symbol, Symbol> ObjectTable;


ClassTree class_tree_;
MethodTable *method_table;
ObjectTable *object_table;
Class_ this_class_;

int errors() { return semant_errors; }
ostream& semant_error();
ostream& semant_error(Class_ c);
ostream& semant_error(Symbol filename, tree_node *t);

// own code
bool CheckAndBuildTypeTree();
void ScopeCheckingAndTypeChecking();

void install_basic_classes() {

    // The tree package uses these globals to annotate the classes built below.
   // curr_lineno  = 0;
    Symbol filename = stringtable.add_string("<basic class>");
    
    // The following demonstrates how to create dummy parse trees to
    // refer to basic Cool classes.  There's no need for method
    // bodies -- these are already built into the runtime system.
    
    // IMPORTANT: The results of the following expressions are
    // stored in local variables.  You will want to do something
    // with those variables at the end of this method to make this
    // code meaningful.

    // 
    // The Object class has no parent class. Its methods are
    //        abort() : Object    aborts the program
    //        type_name() : Str   returns a string representation of class name
    //        copy() : SELF_TYPE  returns a copy of the object
    //
    // There is no need for method bodies in the basic classes---these
    // are already built in to the runtime system.

    Class_ Object_class =
	class_(Object, 
	       No_class,
	       append_Features(
			       append_Features(
					       single_Features(method(cool_abort, nil_Formals(), Object, no_expr())),
					       single_Features(method(type_name, nil_Formals(), Str, no_expr()))),
			       single_Features(method(copy, nil_Formals(), SELF_TYPE, no_expr()))),
	       filename);

    // 
    // The IO class inherits from Object. Its methods are
    //        out_string(Str) : SELF_TYPE       writes a string to the output
    //        out_int(Int) : SELF_TYPE            "    an int    "  "     "
    //        in_string() : Str                 reads a string from the input
    //        in_int() : Int                      "   an int     "  "     "
    //
    Class_ IO_class = 
	class_(IO, 
	       Object,
	       append_Features(
			       append_Features(
					       append_Features(
							       single_Features(method(out_string, single_Formals(formal(arg, Str)),
										      SELF_TYPE, no_expr())),
							       single_Features(method(out_int, single_Formals(formal(arg, Int)),
										      SELF_TYPE, no_expr()))),
					       single_Features(method(in_string, nil_Formals(), Str, no_expr()))),
			       single_Features(method(in_int, nil_Formals(), Int, no_expr()))),
	       filename);  

    //
    // The Int class has no methods and only a single attribute, the
    // "val" for the integer. 
    //
    Class_ Int_class =
	class_(Int, 
	       Object,
	       single_Features(attr(val, prim_slot, no_expr())),
	       filename);

    //
    // Bool also has only the "val" slot.
    //
    Class_ Bool_class =
	class_(Bool, Object, single_Features(attr(val, prim_slot, no_expr())),filename);

    //
    // The class Str has a number of slots and operations:
    //       val                                  the length of the string
    //       str_field                            the string itself
    //       length() : Int                       returns length of the string
    //       concat(arg: Str) : Str               performs string concatenation
    //       substr(arg: Int, arg2: Int): Str     substring selection
    //       
    Class_ Str_class =
	class_(Str, 
	       Object,
	       append_Features(
			       append_Features(
					       append_Features(
							       append_Features(
									       single_Features(attr(val, Int, no_expr())),
									       single_Features(attr(str_field, prim_slot, no_expr()))),
							       single_Features(method(length, nil_Formals(), Int, no_expr()))),
					       single_Features(method(concat, 
								      single_Formals(formal(arg, Str)),
								      Str, 
								      no_expr()))),
			       single_Features(method(substr, 
						      append_Formals(single_Formals(formal(arg, Int)), 
								     single_Formals(formal(arg2, Int))),
						      Str, 
						      no_expr()))),
	       filename);

	// add basic class to classes
	Class_ basic_classes[] = {Object_class, IO_class, Int_class, Bool_class, Str_class};	
	int basic_classes_size = sizeof(basic_classes) / sizeof(Class_);
	for (int i = basic_classes_size - 1; i >= 0; --i) {
		classes_ = append_Classes(single_Classes(basic_classes[i]), classes_);
	}

}


ostream& semant_error(Class_ c)
{                                                             
    return semant_error(c->get_filename(),c);
}    

ostream& semant_error(Symbol filename, tree_node *t)
{
    error_stream << filename << ":" << t->get_line_number() << ": ";
    return semant_error();
}

ostream& semant_error()                  
{                                                 
    semant_errors++;                            
    return error_stream;
} 


/*   This is the entry point to the semantic checker.

     Your checker should do the following two things:

     1) Check that the program is semantically correct
     2) Decorate the abstract syntax tree with type information
        by setting the `type' field in each Expression node.
        (see `tree.h')

     You are free to first do 1), make sure you catch all semantic
     errors. Part 2) can be done in a second stage, when you want
     to build mycoolc.
 */
void program_class::semant()
{
    initialize_constants();

    classes_ = classes;
    method_table = new MethodTable();
	object_table = new ObjectTable();

	install_basic_classes();
    
    /* some semantic analysis code may go here */

	if (CheckAndBuildTypeTree()) {
		ScopeCheckingAndTypeChecking();
	}

    delete object_table;
    delete method_table;

    if (errors()) {
        cerr << "Compilation halted due to static semantic errors." << endl;
        exit(1);
    }

}

bool CheckAndBuildTypeTree()
{
	bool has_main = false;
	// check class redefined
	for (int i = classes_->first(); classes_->more(i); i = classes_->next(i)) {
		Class_ class_ = classes_->nth(i);

		if (class_->get_parent() == Int || class_->get_parent() == Bool || class_->get_parent() == Str) {
			semant_error(class_) << "class " << class_->get_name()->get_string() << " inherit from " << class_->get_parent()->get_string() << " is not allowed" << endl;
			return false;
		}
		ClassTreeNode &node = class_tree_[class_->get_name()];
		if (node.inited) {
			semant_error(class_) << "class " << class_->get_name()->get_string() << " redefined" << endl;
			return false;
		} 

		if (class_->get_name() == Main) {
			has_main = true;
		}

		if (class_->get_name() == SELF_TYPE) {
			semant_error(class_) << "Redefinition of basic class SELF_TYPE." << endl;
			return false;
		}

		node.inited = true;
		node.parent = class_->get_parent();
		node.class_ = class_;

		if (class_->get_parent() != No_class) {
			class_tree_[class_->get_parent()].children.push_back(class_->get_name());
		}


	}

	if (!has_main) {
		semant_error() << "Class Main is not defined." << endl;
		return false;
	}


	// check inherited class defined

	for (int i = classes_->first(); classes_->more(i); i = classes_->next(i)) {
		Class_ class_ = classes_->nth(i);
		if (class_->get_parent() != No_class && !class_tree_[class_->get_parent()].inited) {
			semant_error(class_) << "the parent " << class_->get_parent()->get_string() << " of class " << class_->get_name()->get_string() << " is not defined" << endl;
			return false;
		}
	}


	// check cycle
	std::map<Symbol, int> visit_state;
	for (ClassTree::const_iterator it = class_tree_.begin(); it != class_tree_.end(); ++it) {
		if (0 == visit_state[it->first] && HasCycle(it->first, visit_state)) {
			semant_error(it->second.class_) << "cycle inheritance" << endl;
			return false;
		}
	}

	return true;
}

bool HasCycle(Symbol s, std::map<Symbol, int> &visit_state)
{
	visit_state[s] = -1;

	const std::vector<Symbol> &children = class_tree_[s].children;
	for (int i = 0; i < static_cast<int>(children.size()); i++) {
		Symbol t = children[i];
		if (-1 == visit_state[t]) {
			return true;
		} else if (0 == visit_state[t] && HasCycle(t, visit_state)) {
			return true;
		}
	}
	visit_state[s] = 1;

	return false;
}

int CheckTypeClass(Symbol type)
{
	if (class_tree_.find(type) != class_tree_.end()) {
		return TYPE_CLASS;
	}
   	if (SELF_TYPE == type) {
		return TYPE_SELF;
	}
	if ('_' == type->get_string()[0]) {
		return TYPE_INSIDE;
	}
	return 0;

}

void ScopeCheckingAndTypeChecking()
{
	method_table->enterscope();
	CreateMethodTable(Object, std::vector<Feature>());
	TraversalClass(Object);
	method_table->exitscope();
}





Symbol SetExprType(Expression expr)
{
    Symbol type = expr->ExprType();
    expr->set_type(type);
	return type;

}

Symbol int_const_class::ExprType()
{
    return Int;
}

Symbol bool_const_class::ExprType()
{
    return Bool;
}

Symbol string_const_class::ExprType()
{
    return Str;
}

Symbol object_class::ExprType()
{
	Symbol type = Object, *object_type = object_table->lookup(name);
	if (NULL != object_type) {
		type = *object_type;
	} else {
		semant_error(this_class_) << "identifier " << name->get_string() << " is not defined" << endl;
	}
	return type;
}

Symbol assign_class::ExprType()
{
    Symbol type = SetExprType(expr);
    if (self == name) {
        semant_error(this_class_) << "assign to self is illegal" << endl;
    } else {
        Symbol *object_type = object_table->lookup(name);
        if (NULL != object_type) {
            if (!TypeCheck(*object_type, type)) {
                semant_error(this_class_) << "type " << type->get_string() << " cannot cast to " << (*object_type)->get_string()
                                     << " in assign to identifier " << name->get_string() << endl;
            }
        } else {
            semant_error(this_class_) << "identifier " << name->get_string() << " is not defined" << endl;
        }
    }
    return type;
}

Symbol new__class::ExprType()
{
    Symbol type = Object;
    int ret = CheckTypeClass(type_name);
    if (TYPE_SELF == ret || TYPE_CLASS == ret) {
        type = type_name;
    } else {
        semant_error(this_class_) << "new expr type " << type_name->get_string() << " is not defined" << endl;
    }
    return type;
}

Symbol dispatch_class::ExprType()
{
    Symbol type = SetExprType(expr);
    return CommonDispatchExpr(type, type, name, actual);
}

Symbol static_dispatch_class::ExprType()
{
	return CommonDispatchExpr(SetExprType(expr), type_name, name, actual);
}


Symbol CommonDispatchExpr(Symbol t0, Symbol t1, Symbol name, Expressions actual)
{
	Symbol type = Object;

	std::vector<Symbol> actual_types;
	for (int i = actual->first(); actual->more(i); i = actual->next(i)) {
		actual_types.push_back(SetExprType(actual->nth(i)));
	}

	if (CheckTypeClass(t1)) {
		if (TypeCheck(t1, t0)) {
			MethodTableKey key = std::make_pair(Fmt_SELF_TYPE(t1), name);
			MethodTableValue *val = method_table->lookup(key);
			if (NULL == val) {
				semant_error(this_class_) << "dispatch method " << name->get_string() << " is not defined" << endl;
			} else {
				Formals formals = val->first;
				if (formals->len() == static_cast<int>(actual_types.size())) {
					int i = 0;
					int j = formals->first();
					while (i < static_cast<int>(actual_types.size()) && formals->more(j)) {
						if (!TypeCheck(formals->nth(j)->get_type_decl(), actual_types[i])) {
							semant_error(this_class_) << "the actual param type " << actual_types[i]->get_string() 
                                                 << " of dispatch method " << name->get_string() << " cannot cast to " << formals->nth(j)->get_type_decl()->get_string() << endl;
						}
						i++;
						j = formals->next(j);
					}

				} else {
					semant_error(this_class_) << "length of dispatch method " << name->get_string() 
                                         << " actual param is not same with formal param" << endl;
				}
				if (SELF_TYPE != val->second) {
					type = val->second;
				} else {
					type = t0;
				}
			}
		} else {
			semant_error(this_class_) << "static dispatch method " << name->get_string()
                                 << " in type " << t1->get_string() << " is illegal" << endl;
		}
	} else {
		semant_error(this_class_) << "the type " << t1->get_string() << " in static dispatch method "
                             << name->get_string() << " is illegal"<< endl;
	}

	return type;

}

Symbol cond_class::ExprType()
{
    Symbol pred_type = SetExprType(pred);
	if (Bool != pred_type) {
		semant_error(this_class_) << "cond pred type is wrong" << endl;
	}	
	Symbol then_exp_type = SetExprType(then_exp),
		   else_exp_type = SetExprType(else_exp);
	return Lub(then_exp_type, else_exp_type);
}

Symbol block_class::ExprType()
{
    Symbol type;
	for (int i = body->first(); body->more(i); i = body->next(i)) {
		type = SetExprType(body->nth(i));
	}
	return type;
}

Symbol let_class::ExprType()
{
	Symbol type, type_decl_1 = Object;

	if (CheckTypeClass(type_decl)) {
		type_decl_1 = type_decl;
	} else {
		semant_error(this_class_) << "let expr the type " << type_decl->get_string() << " is not defined" << endl;
	}
	Symbol init_type = SetExprType(init);
	if (!TypeCheck(type_decl_1, init_type)) {
		semant_error(this_class_) << "let expr the init type " << init_type->get_string()
                             << " cannot cast to " << type_decl_1->get_string() << endl;
	}
	object_table->enterscope();
	if (self != identifier) {  
		object_table->addid(identifier, new Symbol(type_decl_1));
	} else {
		semant_error(this_class_) << "let expr bind self is illegal" << endl;
	}
	type = SetExprType(body);
	object_table->exitscope();
	return type;
}


Symbol typcase_class::ExprType()
{
	Symbol type = Object;
	bool inited = false;
	SetExprType(expr);

	std::set<Symbol> type_decls;
	for (int i = cases->first(); cases->more(i); i = cases->next(i)) {
		branch_class *branch = static_cast<branch_class*>(cases->nth(i));
		object_table->enterscope();
		Symbol type_decl = Object;
		if (TYPE_CLASS == CheckTypeClass(branch->get_type_decl())) {
			type_decl = branch->get_type_decl();
			if (type_decls.find(type_decl) == type_decls.end()) {
				type_decls.insert(type_decl);
			} else {
				semant_error(this_class_) << "case expr the type decl " << branch->get_type_decl()->get_string() 
                                     << " is conflict" << endl; 
			}
		} else {
			semant_error(this_class_) << "case expr the type " << branch->get_type_decl()->get_string() 
                                 << " is not defined" << endl; 
		}
		object_table->addid(branch->get_name(), new Symbol(type_decl));

		if (!inited)  {
			inited = true;
			type = SetExprType(branch->get_expr());
		} else {
			type = Lub(type, SetExprType(branch->get_expr())); 
		}

		object_table->exitscope();
		
	}
	return type;
}

Symbol loop_class::ExprType()
{
    Symbol type = Object;
	if (SetExprType(pred) != Bool) {
		semant_error(this_class_) << "the type of pred in while is not bool" << endl;
	}
	SetExprType(body);
	return type;
}

Symbol isvoid_class::ExprType()
{
    Symbol type = Bool;

	SetExprType(e1);

	return type;
}

bool InlineClass(Symbol type)
{
    return type == Int || type == Str || type == Bool;
}

Symbol eq_class::ExprType()
{
	Symbol type = Bool;
	Symbol e1_type = SetExprType(e1), e2_type = SetExprType(e2);
	if (InlineClass(e1_type) || InlineClass(e2_type)) {
		if (e1_type != e2_type) {
			semant_error(this_class_) << "eq expr type of e1 must be same with type of e2, when one type is Int, String, Bool" << endl;
		}
	}
	return type;
}

Symbol UnaryExpr(Expression e1, Symbol target_type)
{
	Symbol type = target_type;
	if (SetExprType(e1) != target_type) {
		semant_error(this_class_) << "the type of unary expr is not " <<  target_type->get_string() << endl;
	}
	return type;
}


Symbol BinaryExpr(Expression e1, Expression e2, Symbol target_type)
{
	Symbol type = target_type;

	Symbol e1_type = SetExprType(e1), e2_type = SetExprType(e2);

	if (e1_type != Int) {
		semant_error(this_class_) << "the type of binary expr e1 is not allowed" << endl;
	} 
	if (e2_type != Int) {
		semant_error(this_class_) << "the type of binary expr e2 is not allowed" << endl;
	} 
	return type;
}



Symbol Fmt_SELF_TYPE(Symbol type)
{
    if (type == SELF_TYPE) {
        return this_class_->get_name();
    } 
    return type;
}

Symbol comp_class::ExprType()
{
    return UnaryExpr(e1, Bool);
}

Symbol neg_class::ExprType()
{
    return UnaryExpr(e1, Int);
}

Symbol lt_class::ExprType()
{
    return BinaryExpr(e1, e2, Bool);
}

Symbol leq_class::ExprType()
{
    return BinaryExpr(e1, e2, Bool);
}

Symbol plus_class::ExprType()
{
    return BinaryExpr(e1, e2, Int);
}

Symbol sub_class::ExprType()
{
    return BinaryExpr(e1, e2, Int);
}

Symbol mul_class::ExprType()
{
    return BinaryExpr(e1, e2, Int);
}

Symbol divide_class::ExprType()
{
    return BinaryExpr(e1, e2, Int);
}

Symbol no_expr_class::ExprType()
{
    return No_type;
}

Symbol Lub(Symbol t1, Symbol t2)
{
	bool least_r1 = false, least_r2 = false;
	t1 = RegularType(t1, least_r1);
	t2 = RegularType(t2, least_r2);
	if (least_r1) {
		return t2;
	} else {
		if (least_r2) {
			return t1;
		} else {
			if (t1 == SELF_TYPE && t2 == SELF_TYPE) {
				return SELF_TYPE;
			} else {
				if (t1 == SELF_TYPE) {
					t1 = this_class_->get_name();
				}

				if (t2 == SELF_TYPE) {
					t2 = this_class_->get_name();
				}
				Symbol s1 = t1;
				while (s1 != No_class) {
					Symbol s2 = t2;
					while (s2 != No_class) {
						if (s1 == s2) {
							return s1;
						}
						s2 = class_tree_[s2].parent;
					}
					s1 = class_tree_[s1].parent;
				}
			}
		}
	}
	// never reach this
	return Object;
}


bool TypeCheck(Symbol t1, Symbol t2)
{
	bool least_r1 = false, least_r2 = false;
	t1 = RegularType(t1, least_r1);
	t2 = RegularType(t2, least_r2);
	if (least_r1) {
		return least_r2;
	} else {
		if (least_r2) {
			return true;
		} else {
			if (t1 == SELF_TYPE) {
				if (t2 == SELF_TYPE) {
					return true;
				} else {
					return false;
				}
			} else {
				if (t2 == SELF_TYPE) {
					t2 = this_class_->get_name();
				}
				Symbol s = t2;
				while (s != No_class) {
					if (s == t1) {
						return true;
					}
					s = class_tree_[s].parent;
				}
			}

		}
	}
	return false;
}

Symbol RegularType(Symbol t, bool &least)
{
	int ret = CheckTypeClass(t);
	if (TYPE_CLASS == ret || TYPE_SELF == ret) {
        return t;
	} else if (TYPE_INSIDE == ret) {
		least = true;
		return t;
	} else {
        return Object;
	}
}

void CreateMethodTable(Symbol curr_class, const std::vector<Feature> &ancestor_method)
{
	const ClassTreeNode &node = class_tree_[curr_class];
    this_class_ = node.class_;
	std::vector<Feature> new_method;

	Features features = this_class_->get_features();
	for (int i = features->first(); features->more(i); i = features->next(i)) {
		Feature feature = features->nth(i);
		if (feature->GetType() == Feature_class::METHOD) {
            if (feature->CheckTypeAndInsertToSymbolTab()) {
                new_method.push_back(feature);
            }
		}
	}

	std::vector<Feature> next_ancestor_method;

	for (int i = 0; i < static_cast<int>(ancestor_method.size()); i++) {
        if (ancestor_method[i]->AncestorAdd()) {
            next_ancestor_method.push_back(ancestor_method[i]);
        }
	}
    
	for (int i = 0; i < static_cast<int>(new_method.size()); i++) {
		next_ancestor_method.push_back(new_method[i]);
	}

	for (int i = 0; i < static_cast<int>(node.children.size()); i++) {
		CreateMethodTable(node.children[i], next_ancestor_method);
	}
}

void TraversalClass(Symbol curr_class) 
{
    const ClassTreeNode &node = class_tree_[curr_class];
	this_class_ = node.class_;

	object_table->enterscope();

	Features features = this_class_->get_features();
	// add attr to object_table and check type
	for (int i = features->first(); features->more(i); i = features->next(i)) {
		Feature feature = features->nth(i);
		if (feature->GetType() == Feature_class::ATTR) {
            feature->CheckTypeAndInsertToSymbolTab();

            feature->AncestorAdd();
            

		}
	}

	object_table->addid(self, new Symbol(SELF_TYPE));

	// type check attr expression
	for (int i = features->first(); features->more(i); i = features->next(i)) {
		Feature feature = features->nth(i);
		if (feature->GetType() == Feature_class::ATTR) {
            feature->CheckExprType();
		}
	}

    // i think check and output attr first, and then check method is more friendly
	for (int i = features->first(); features->more(i); i = features->next(i)) {
		Feature feature = features->nth(i);
		if (feature->GetType() == Feature_class::METHOD) {
            feature->CheckExprType();
		}
	}


	for (int i = 0; i < static_cast<int>(node.children.size()); i++) {
		TraversalClass(node.children[i]);
	}


	object_table->exitscope();
}

bool attr_class::CheckTypeAndInsertToSymbolTab()
{
    Symbol attr_type = Object;
    if (CheckTypeClass(type_decl)) {
        attr_type = type_decl;
    } else {
        semant_error(this_class_) << "the type " << type_decl->get_string() << " of attribute "
                             << name->get_string() << " is not defined" << endl;
    }
    if (name != self) {
        if (NULL == object_table->lookup(name)) {
            object_table->addid(name, new Symbol(attr_type));
            return true;
        } else {
            if (NULL != object_table->probe(name)) {
                semant_error(this_class_) << "the attribute " << name->get_string() << " is redefined" << endl;
            } else {
                semant_error(this_class_) << "the attribute " << name->get_string()
                                          <<  " is inherited attributes, cannot redefined" << endl;
            }
        }
    } else {
        semant_error(this_class_) << "the attribute cannot be self" << endl;
    }
    return false;

}

bool attr_class::AncestorAdd()
{
    // implicit add in object_table
    return true;
}

void attr_class::CheckExprType()
{
    if (!TypeCheck(type_decl, SetExprType(init))) {
        semant_error(this_class_) << "the init of the attribute " << name->get_string()
                             << " is not cast to " <<  type_decl->get_string() << endl;
    }
}

bool method_class::CheckTypeAndInsertToSymbolTab()
{
    MethodTableKey key = std::make_pair(this_class_->get_name(), name);
    if (method_table->probe(key) == NULL) {
        if (!CheckTypeClass(return_type)) {
            semant_error(this_class_) << "return type " << return_type->get_string() << " of method "
                                 << name->get_string() << " is not defined" << endl;
        }
        for (int j = formals->first(); formals->more(j); j = formals->next(j)) {
            Formal formal = formals->nth(j);
            int ret = CheckTypeClass(formal->get_type_decl());
            if (0 == ret || TYPE_SELF == ret) {
                std::string msg;
                if (0 == ret) {
                    msg = " is not defined";
                } else if (TYPE_SELF == ret) {
                    msg = " is SELF_TYPE, not allowed";
                }
                semant_error(this_class_) << "formal parameter type " << formal->get_type_decl()->get_string() << " of method " 
                                          << name->get_string() << msg << endl;
            } 				
        }
        method_table->addid(key, new MethodTableValue(formals, return_type));
        return true;
    } else {
        semant_error(this_class_) << "method " << name->get_string() << " is redefined" << endl;
        return false;
    }
}

bool method_class::AncestorAdd()
{
    MethodTableKey key = std::make_pair(this_class_->get_name(), name);
    MethodTableValue *val = method_table->probe(key);
    if (NULL == val) {
        method_table->addid(key, new MethodTableValue(formals, return_type));
        return true;
    } else {
        int is_valid = 0;
        if (return_type == val->second) {
            Formals now_formals = val->first;
            if (formals->len() == now_formals->len()) {
                int j = formals->first();
                int k = now_formals->first();
                while (formals->more(j) && now_formals->more(k)) {
                    if (formals->nth(j)->get_type_decl() != now_formals->nth(k)->get_type_decl()) {
                        is_valid = -1;
                        break;
                    }
                    j = formals->next(j);
                    k = now_formals->next(k);
                }

            } else {
                is_valid = -2;
            }
        } else {
            is_valid = -3;
        }
        if (0 != is_valid) {
            std::string msg;
            if (-1 == is_valid) {
                msg = "the type of formal parameter is not same";
            } else if (-2 == is_valid) {
                msg = "the number of formal parameters is not same";
            } else if (-3 == is_valid) {
                msg = "return type is not same";
            }
            semant_error(this_class_) << "the redefinition of method " << name->get_string() << " in class " << this_class_->get_name()->get_string() << " is invalid: " << msg << endl;
        }
        return false;
    }    
    
}

void method_class::CheckExprType()
{
    object_table->enterscope();

    for (int j = formals->first(); formals->more(j); j = formals->next(j)) {
        Formal formal = formals->nth(j);
        if (NULL == object_table->probe(formal->get_name())) {
            Symbol formal_type = Object;
            int ret = CheckTypeClass(formal->get_type_decl());
            if (0 != ret && TYPE_SELF != ret) {
                formal_type = formal->get_type_decl();
            }
            if (formal->get_name() != self) {
                object_table->addid(formal->get_name(), new Symbol(formal_type));
            } else {
                semant_error(this_class_) << "formal param of method "
                                          << name->get_string() << " cannot be self" << endl;
            }
        } else {
            semant_error(this_class_) << "the formal param " << formal->get_name()->get_string() << " of the method " 
                                 << name->get_string() << " is redefined" << endl;
        }

    }

    if (!TypeCheck(return_type, SetExprType(expr))) {
        semant_error(this_class_) << "the expr type of the method " << name->get_string() << " is not cast to " << return_type->get_string() << endl;
    }


    object_table->exitscope();
}







