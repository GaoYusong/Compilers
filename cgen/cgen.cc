
//**************************************************************
//
// Code generator SKELETON
//
// Read the comments carefully. Make sure to
//    initialize the base class tags in
//       `CgenClassTable::CgenClassTable'
//
//    Add the label for the dispatch tables to
//       `IntEntry::code_def'
//       `StringEntry::code_def'
//       `BoolConst::code_def'
//
//    Add code to emit everyting else that is needed
//       in `CgenClassTable::code'
//
//
// The files as provided will produce code to begin the code
// segments, declare globals, and emit constants.  You must
// fill in the rest.
//
//**************************************************************

#include "cgen.h"
#include "cgen_gc.h"
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <iostream>
#include <sstream>
#include <assert.h>



extern void emit_string_constant(ostream& str, char *s);
extern int cgen_debug;

//
// Three symbols from the semantic analyzer (semant.cc) are used.
// If e : No_type, then no code is generated for e.
// Special code is generated for new SELF_TYPE.
// The name "self" also generates code different from other references.
//
//////////////////////////////////////////////////////////////////////
//
// Symbols
//
// For convenience, a large number of symbols are predefined here.
// These symbols include the primitive type and method names, as well
// as fixed names used by the runtime system.
//
//////////////////////////////////////////////////////////////////////
Symbol 
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

static char *gc_init_names[] =
  { "_NoGC_Init", "_GenGC_Init", "_ScnGC_Init" };
static char *gc_collect_names[] =
  { "_NoGC_Collect", "_GenGC_Collect", "_ScnGC_Collect" };


//  BoolConst is a class that implements code generation for operations
//  on the two booleans, which are given global names here.
BoolConst falsebool(FALSE);
BoolConst truebool(TRUE);

SymbolTable<Symbol, int> *locations = NULL;
CgenNodeP self_ = NULL;
CgenClassTableP class_table = NULL;
int sp_level;
int LABEL_COUNT = 0;
std::map<Symbol, int> n2tag;

#define CLASSINHERITTAB "class_inheritTab"

//*********************************************************
//
// Define method for code generation
//
// This is the method called by the compiler driver
// `cgtest.cc'. cgen takes an `ostream' to which the assembly will be
// emmitted, and it passes this and the class list of the
// code generator tree to the constructor for `CgenClassTable'.
// That constructor performs all of the work of the code
// generator.
//
//*********************************************************

void program_class::cgen(ostream &os) 
{
  // spim wants comments to start with '#'
  os << "# start of generated code\n";

  initialize_constants();
  CgenClassTable *codegen_classtable = new CgenClassTable(classes,os);

  os << "\n# end of generated code\n";
}


//////////////////////////////////////////////////////////////////////////////
//
//  emit_* procedures
//
//  emit_X  writes code for operation "X" to the output stream.
//  There is an emit_X for each opcode X, as well as emit_ functions
//  for generating names according to the naming conventions (see emit.h)
//  and calls to support functions defined in the trap handler.
//
//  Register names and addresses are passed as strings.  See `emit.h'
//  for symbolic names you can use to refer to the strings.
//
//////////////////////////////////////////////////////////////////////////////

static void emit_load(const char *dest_reg, int offset, const char *source_reg, ostream& s)
{
  s << LW << dest_reg << " " << offset * WORD_SIZE << "(" << source_reg << ")" 
    << endl;
}

static void emit_store(const char *source_reg, int offset, const char *dest_reg, ostream& s)
{
  s << SW << source_reg << " " << offset * WORD_SIZE << "(" << dest_reg << ")"
      << endl;
}

static void emit_load_imm(char *dest_reg, int val, ostream& s)
{ s << LI << dest_reg << " " << val << endl; }

static void emit_load_address(const char *dest_reg, const char *address, ostream& s)
{ s << LA << dest_reg << " " << address << endl; }

static void emit_partial_load_address(char *dest_reg, ostream& s)
{ s << LA << dest_reg << " "; }

static void emit_load_bool(char *dest, const BoolConst& b, ostream& s)
{
  emit_partial_load_address(dest,s);
  b.code_ref(s);
  s << endl;
}

static void emit_load_string(char *dest, StringEntry *str, ostream& s)
{
  emit_partial_load_address(dest,s);
  str->code_ref(s);
  s << endl;
}

static void emit_load_int(char *dest, IntEntry *i, ostream& s)
{
  emit_partial_load_address(dest,s);
  i->code_ref(s);
  s << endl;
}

static void emit_move(const char *dest_reg, const char *source_reg, ostream& s)
{ s << MOVE << dest_reg << " " << source_reg << endl; }

static void emit_neg(char *dest, char *src1, ostream& s)
{ s << NEG << dest << " " << src1 << endl; }

static void emit_add(char *dest, char *src1, char *src2, ostream& s)
{ s << ADD << dest << " " << src1 << " " << src2 << endl; }

static void emit_addu(char *dest, char *src1, char *src2, ostream& s)
{ s << ADDU << dest << " " << src1 << " " << src2 << endl; }

static void emit_addiu(const char *dest, const char *src1, int imm, ostream& s)
{ s << ADDIU << dest << " " << src1 << " " << imm << endl; }

static void emit_div(char *dest, char *src1, char *src2, ostream& s)
{ s << DIV << dest << " " << src1 << " " << src2 << endl; }

static void emit_mul(char *dest, char *src1, char *src2, ostream& s)
{ s << MUL << dest << " " << src1 << " " << src2 << endl; }

static void emit_sub(char *dest, char *src1, char *src2, ostream& s)
{ s << SUB << dest << " " << src1 << " " << src2 << endl; }

static void emit_sll(char *dest, char *src1, int num, ostream& s)
{ s << SLL << dest << " " << src1 << " " << num << endl; }

static void emit_jalr(char *dest, ostream& s)
{ s << JALR << "\t" << dest << endl; }

static void emit_jal(const char *address,ostream &s)
{ s << JAL << address << endl; }

static void emit_return(ostream& s)
{ s << RET << endl; }

static void emit_gc_assign(ostream& s)
{ s << JAL << "_GenGC_Assign" << endl; }

static void emit_disptable_ref(Symbol sym, ostream& s)
{  s << sym << DISPTAB_SUFFIX; }

static void emit_init_ref(Symbol sym, ostream& s)
{ s << sym << CLASSINIT_SUFFIX; }

static void emit_label_ref(int l, ostream &s)
{ s << "label" << l; }

static void emit_protobj_ref(Symbol sym, ostream& s)
{ s << sym << PROTOBJ_SUFFIX; }

static void emit_method_ref(Symbol classname, Symbol methodname, ostream& s)
{ s << classname << METHOD_SEP << methodname; }

static void emit_label_def(int l, ostream &s)
{
  emit_label_ref(l,s);
  s << ":" << endl;
}

static void emit_beqz(char *source, int label, ostream &s)
{
  s << BEQZ << source << " ";
  emit_label_ref(label,s);
  s << endl;
}

static void emit_beq(char *src1, char *src2, int label, ostream &s)
{
  s << BEQ << src1 << " " << src2 << " ";
  emit_label_ref(label,s);
  s << endl;
}

static void emit_bne(char *src1, char *src2, int label, ostream &s)
{
  s << BNE << src1 << " " << src2 << " ";
  emit_label_ref(label,s);
  s << endl;
}

static void emit_bleq(char *src1, char *src2, int label, ostream &s)
{
  s << BLEQ << src1 << " " << src2 << " ";
  emit_label_ref(label,s);
  s << endl;
}

static void emit_blt(char *src1, char *src2, int label, ostream &s)
{
  s << BLT << src1 << " " << src2 << " ";
  emit_label_ref(label,s);
  s << endl;
}

static void emit_blti(char *src1, int imm, int label, ostream &s)
{
  s << BLT << src1 << " " << imm << " ";
  emit_label_ref(label,s);
  s << endl;
}

static void emit_bgti(char *src1, int imm, int label, ostream &s)
{
  s << BGT << src1 << " " << imm << " ";
  emit_label_ref(label,s);
  s << endl;
}

static void emit_branch(int l, ostream& s)
{
  s << BRANCH;
  emit_label_ref(l,s);
  s << endl;
}

//
// Push a register on the stack. The stack grows towards smaller addresses.
//
static void emit_push(char *reg, ostream& str)
{
  emit_store(reg,0,SP,str);
  emit_addiu(SP,SP,-4,str);
}

//
// Fetch the integer value in an Int object.
// Emits code to fetch the integer value of the Integer object pointed
// to by register source into the register dest
//
static void emit_fetch_int(char *dest, char *source, ostream& s)
{ emit_load(dest, DEFAULT_OBJFIELDS, source, s); }

//
// Emits code to store the integer value contained in register source
// into the Integer object pointed to by dest.
//
static void emit_store_int(char *source, char *dest, ostream& s)
{ emit_store(source, DEFAULT_OBJFIELDS, dest, s); }


static void emit_test_collector(ostream &s)
{
  emit_push(ACC, s);
  emit_move(ACC, SP, s); // stack end
  emit_move(A1, ZERO, s); // allocate nothing
  s << JAL << gc_collect_names[cgen_Memmgr] << endl;
  emit_addiu(SP,SP,4,s);
  emit_load(ACC,0,SP,s);
}

static void emit_gc_check(char *source, ostream &s)
{
  if (source != (char*)A1) emit_move(A1, source, s);
  s << JAL << "_gc_check" << endl;
}


///////////////////////////////////////////////////////////////////////////////
//
// coding strings, ints, and booleans
//
// Cool has three kinds of constants: strings, ints, and booleans.
// This section defines code generation for each type.
//
// All string constants are listed in the global "stringtable" and have
// type StringEntry.  StringEntry methods are defined both for String
// constant definitions and references.
//
// All integer constants are listed in the global "inttable" and have
// type IntEntry.  IntEntry methods are defined for Int
// constant definitions and references.
//
// Since there are only two Bool values, there is no need for a table.
// The two booleans are represented by instances of the class BoolConst,
// which defines the definition and reference methods for Bools.
//
///////////////////////////////////////////////////////////////////////////////

//
// Strings
//
void StringEntry::code_ref(ostream& s)
{
  s << STRCONST_PREFIX << index;
}

//
// Emit code for a constant String.
// You should fill in the code naming the dispatch table.
//

void StringEntry::code_def(ostream& s, int stringclasstag)
{
  IntEntryP lensym = inttable.add_int(len);

  // Add -1 eye catcher
  s << WORD << "-1" << endl;

  code_ref(s);  s  << LABEL                                             // label
      << WORD << stringclasstag << endl                                 // tag
      << WORD << (DEFAULT_OBJFIELDS + STRING_SLOTS + (len+4)/4) << endl // size
      << WORD;


 /***** Add dispatch information for class String ******/

      emit_disptable_ref(Str, s); s << endl;       // dispatch table
      s << WORD;  lensym->code_ref(s);  s << endl;            // string length
  emit_string_constant(s,str);                                // ascii string
  s << ALIGN;                                                 // align to word
}

//
// StrTable::code_string
// Generate a string object definition for every string constant in the 
// stringtable.
//
void StrTable::code_string_table(ostream& s, int stringclasstag)
{  
  for (List<StringEntry> *l = tbl; l; l = l->tl())
    l->hd()->code_def(s,stringclasstag);
}

//
// Ints
//
void IntEntry::code_ref(ostream &s)
{
  s << INTCONST_PREFIX << index;
}

//
// Emit code for a constant Integer.
// You should fill in the code naming the dispatch table.
//

void IntEntry::code_def(ostream &s, int intclasstag)
{
  // Add -1 eye catcher
  s << WORD << "-1" << endl;

  code_ref(s);  s << LABEL                                // label
      << WORD << intclasstag << endl                      // class tag
      << WORD << (DEFAULT_OBJFIELDS + INT_SLOTS) << endl  // object size
      << WORD; 

 /***** Add dispatch information for class Int ******/
 
      emit_disptable_ref(Int, s); s << endl;    // dispatch table
      s << WORD << str << endl;                           // integer value
}


//
// IntTable::code_string_table
// Generate an Int object definition for every Int constant in the
// inttable.
//
void IntTable::code_string_table(ostream &s, int intclasstag)
{
  for (List<IntEntry> *l = tbl; l; l = l->tl())
    l->hd()->code_def(s,intclasstag);
}


//
// Bools
//
BoolConst::BoolConst(int i) : val(i) { assert(i == 0 || i == 1); }

void BoolConst::code_ref(ostream& s) const
{
  s << BOOLCONST_PREFIX << val;
}
  
//
// Emit code for a constant Bool.
// You should fill in the code naming the dispatch table.
//

void BoolConst::code_def(ostream& s, int boolclasstag)
{
  // Add -1 eye catcher
  s << WORD << "-1" << endl;

  code_ref(s);  s << LABEL                                  // label
      << WORD << boolclasstag << endl                       // class tag
      << WORD << (DEFAULT_OBJFIELDS + BOOL_SLOTS) << endl   // object size
      << WORD;

 /***** Add dispatch information for class Bool ******/
 
      emit_disptable_ref(Bool, s); s << endl;    // dispatch table
      s << WORD << val << endl;                             // value (0 or 1)
}

//////////////////////////////////////////////////////////////////////////////
//
//  CgenClassTable methods
//
//////////////////////////////////////////////////////////////////////////////


void CgenClassTable::code_class_nameTab()
{
	str << CLASSNAMETAB << LABEL;
	for (List<CgenNode> *l = nds; l; l = l->tl())
		l->hd()->code_classname(str);
}

void CgenClassTable::code_class_objTab()
{
	str << CLASSOBJTAB << LABEL;
	for (List<CgenNode> *l = nds; l; l = l->tl())
		l->hd()->code_classobj(str);
}

void CgenClassTable::code_class_inheritTab()
{
	str << CLASSINHERITTAB << LABEL;
	int tagid = 0;
	for (List<CgenNode> *l = nds; l; l = l->tl()) {
		n2tag[l->hd()->name] = tagid++;
	}

	n2tag[No_class] = tagid;


	for (List<CgenNode> *l = nds; l; l = l->tl()) {
		str << WORD << n2tag[l->hd()->parent] << endl;
	}

}

void CgenClassTable::code_class_dispTab()
{
	for (List<CgenNode> *l = nds; l; l = l->tl()) 
		l->hd()->code_dispTab(str);
}

void CgenClassTable::code_class_protObj()
{
	int idx = 0;
	for (List<CgenNode> *l = nds; l; l = l->tl())
		l->hd()->code_protObj(str, idx++);
}

void CgenClassTable::code_class_initObj()
{
	for (List<CgenNode> *l = nds; l; l = l->tl()) 
		l->hd()->code_initObj(str);
}

void CgenClassTable::code_class_method()
{
	for (List<CgenNode> *l = nds; l; l = l->tl())
		l->hd()->code_method(str);
}

//***************************************************
//
//  Emit code to start the .data segment and to
//  declare the global names.
//
//***************************************************

void CgenClassTable::code_global_data()
{
  Symbol main    = idtable.lookup_string(MAINNAME);
  Symbol string  = idtable.lookup_string(STRINGNAME);
  Symbol integer = idtable.lookup_string(INTNAME);
  Symbol boolc   = idtable.lookup_string(BOOLNAME);

  str << "\t.data\n" << ALIGN;
  //
  // The following global names must be defined first.
  //
  str << GLOBAL << CLASSNAMETAB << endl;
  str << GLOBAL; emit_protobj_ref(main,str);    str << endl;
  str << GLOBAL; emit_protobj_ref(integer,str); str << endl;
  str << GLOBAL; emit_protobj_ref(string,str);  str << endl;
  str << GLOBAL; falsebool.code_ref(str);  str << endl;
  str << GLOBAL; truebool.code_ref(str);   str << endl;
  str << GLOBAL << INTTAG << endl;
  str << GLOBAL << BOOLTAG << endl;
  str << GLOBAL << STRINGTAG << endl;

  //
  // We also need to know the tag of the Int, String, and Bool classes
  // during code generation.
  //
  str << INTTAG << LABEL
      << WORD << intclasstag << endl;
  str << BOOLTAG << LABEL 
      << WORD << boolclasstag << endl;
  str << STRINGTAG << LABEL 
      << WORD << stringclasstag << endl;    
}


//***************************************************
//
//  Emit code to start the .text segment and to
//  declare the global names.
//
//***************************************************

void CgenClassTable::code_global_text()
{
  str << GLOBAL << HEAP_START << endl
      << HEAP_START << LABEL 
      << WORD << 0 << endl
      << "\t.text" << endl
      << GLOBAL;
  emit_init_ref(idtable.add_string("Main"), str);
  str << endl << GLOBAL;
  emit_init_ref(idtable.add_string("Int"),str);
  str << endl << GLOBAL;
  emit_init_ref(idtable.add_string("String"),str);
  str << endl << GLOBAL;
  emit_init_ref(idtable.add_string("Bool"),str);
  str << endl << GLOBAL;
  emit_method_ref(idtable.add_string("Main"), idtable.add_string("main"), str);
  str << endl;
}

void CgenClassTable::code_bools(int boolclasstag)
{
  falsebool.code_def(str,boolclasstag);
  truebool.code_def(str,boolclasstag);
}

void CgenClassTable::code_select_gc()
{
  //
  // Generate GC choice constants (pointers to GC functions)
  //
  str << GLOBAL << "_MemMgr_INITIALIZER" << endl;
  str << "_MemMgr_INITIALIZER:" << endl;
  str << WORD << gc_init_names[cgen_Memmgr] << endl;
  str << GLOBAL << "_MemMgr_COLLECTOR" << endl;
  str << "_MemMgr_COLLECTOR:" << endl;
  str << WORD << gc_collect_names[cgen_Memmgr] << endl;
  str << GLOBAL << "_MemMgr_TEST" << endl;
  str << "_MemMgr_TEST:" << endl;
  str << WORD << (cgen_Memmgr_Test == GC_TEST) << endl;
}


//********************************************************
//
// Emit code to reserve space for and initialize all of
// the constants.  Class names should have been added to
// the string table (in the supplied code, is is done
// during the construction of the inheritance graph), and
// code for emitting string constants as a side effect adds
// the string's length to the integer table.  The constants
// are emmitted by running through the stringtable and inttable
// and producing code for each entry.
//
//********************************************************

void CgenClassTable::code_constants()
{
  //
  // Add constants that are required by the code generator.
  //
  stringtable.add_string("");
  inttable.add_string("0");

  stringtable.code_string_table(str,stringclasstag);
  inttable.code_string_table(str,intclasstag);
  code_bools(boolclasstag);
}


CgenClassTable::CgenClassTable(Classes classes, ostream& s) : nds(NULL) , str(s)
{

   enterscope();
   if (cgen_debug) cout << "Building CgenClassTable" << endl;
   install_basic_classes();
   install_classes(classes);
   build_inheritance_tree();

   int idx = 0;
   for (List<CgenNode> *l = nds; l; l = l->tl()) {
	   Symbol name = l->hd()->name;
	   if (name == Str) {
		   stringclasstag = idx;
	   } else if (name == Int) {
		   intclasstag = idx;
	   } else if (name == Bool) {
		   boolclasstag = idx;

	   }
	   idx++;
   }


   locations = new SymbolTable<Symbol, int>();

   class_table = this;

   code();

   delete locations;
   exitscope();
}

void CgenClassTable::install_basic_classes()
{

// The tree package uses these globals to annotate the classes built below.
  //curr_lineno  = 0;
  Symbol filename = stringtable.add_string("<basic class>");

//
// A few special class names are installed in the lookup table but not
// the class list.  Thus, these classes exist, but are not part of the
// inheritance hierarchy.
// No_class serves as the parent of Object and the other special classes.
// SELF_TYPE is the self class; it cannot be redefined or inherited.
// prim_slot is a class known to the code generator.
//
  addid(No_class,
	new CgenNode(class_(No_class,No_class,nil_Features(),filename),
			    Basic,this));
  addid(SELF_TYPE,
	new CgenNode(class_(SELF_TYPE,No_class,nil_Features(),filename),
			    Basic,this));
  addid(prim_slot,
	new CgenNode(class_(prim_slot,No_class,nil_Features(),filename),
			    Basic,this));

// 
// The Object class has no parent class. Its methods are
//        cool_abort() : Object    aborts the program
//        type_name() : Str        returns a string representation of class name
//        copy() : SELF_TYPE       returns a copy of the object
//
// There is no need for method bodies in the basic classes---these
// are already built in to the runtime system.
//
  install_class(
   new CgenNode(
    class_(Object, 
	   No_class,
	   append_Features(
           append_Features(
           single_Features(method(cool_abort, nil_Formals(), Object, no_expr())),
           single_Features(method(type_name, nil_Formals(), Str, no_expr()))),
           single_Features(method(copy, nil_Formals(), SELF_TYPE, no_expr()))),
	   filename),
    Basic,this));

// 
// The IO class inherits from Object. Its methods are
//        out_string(Str) : SELF_TYPE          writes a string to the output
//        out_int(Int) : SELF_TYPE               "    an int    "  "     "
//        in_string() : Str                    reads a string from the input
//        in_int() : Int                         "   an int     "  "     "
//
   install_class(
    new CgenNode(
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
	   filename),	    
    Basic,this));

//
// The Int class has no methods and only a single attribute, the
// "val" for the integer. 
//
   install_class(
    new CgenNode(
     class_(Int, 
	    Object,
            single_Features(attr(val, prim_slot, no_expr())),
	    filename),
     Basic,this));

//
// Bool also has only the "val" slot.
//
    install_class(
     new CgenNode(
      class_(Bool, Object, single_Features(attr(val, prim_slot, no_expr())),filename),
      Basic,this));

//
// The class Str has a number of slots and operations:
//       val                                  ???
//       str_field                            the string itself
//       length() : Int                       length of the string
//       concat(arg: Str) : Str               string concatenation
//       substr(arg: Int, arg2: Int): Str     substring
//       
   install_class(
    new CgenNode(
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
	     filename),
        Basic,this));

}

// CgenClassTable::install_class
// CgenClassTable::install_classes
//
// install_classes enters a list of classes in the symbol table.
//
void CgenClassTable::install_class(CgenNodeP nd)
{
  Symbol name = nd->get_name();

  if (probe(name))
    {
      return;
    }

  // The class name is legal, so add it to the list of classes
  // and the symbol table.
  nds = new List<CgenNode>(nd,nds);
  addid(name,nd);
}

void CgenClassTable::install_classes(Classes cs)
{
  for(int i = cs->first(); cs->more(i); i = cs->next(i))
    install_class(new CgenNode(cs->nth(i),NotBasic,this));
}

//
// CgenClassTable::build_inheritance_tree
//
void CgenClassTable::build_inheritance_tree()
{
  for(List<CgenNode> *l = nds; l; l = l->tl())
      set_relations(l->hd());
}

//
// CgenClassTable::set_relations
//
// Takes a CgenNode and locates its, and its parent's, inheritance nodes
// via the class table.  Parent and child pointers are added as appropriate.
//
void CgenClassTable::set_relations(CgenNodeP nd)
{
  CgenNode *parent_node = probe(nd->get_parent());
  nd->set_parentnd(parent_node);
  parent_node->add_child(nd);
}

void CgenNode::add_child(CgenNodeP n)
{
  children = new List<CgenNode>(n,children);
}

void CgenNode::set_parentnd(CgenNodeP p)
{
  assert(parentnd == NULL);
  assert(p != NULL);
  parentnd = p;
}




void CgenClassTable::code()
{
  if (cgen_debug) cout << "coding global data" << endl;
  code_global_data();

  if (cgen_debug) cout << "choosing gc" << endl;
  code_select_gc();

  if (cgen_debug) cout << "coding constants" << endl;
  code_constants();

//                 Add your code to emit
//                   - prototype objects
//                   - class_nameTab
//                   - dispatch tables
//

  if (cgen_debug) cout << "coding class_nameTab" << endl;
  code_class_nameTab();

  if (cgen_debug) cout << "coding class_objTab" << endl;
  code_class_objTab();

  if (cgen_debug) cout << "coding class_inheritTab" << endl;
  code_class_inheritTab();

  if (cgen_debug) cout << "coding class_dispTab" << endl;
  code_class_dispTab();

  if (cgen_debug) cout << "coding class_protObj" << endl;
  code_class_protObj();



  if (cgen_debug) cout << "coding global text" << endl;
  code_global_text();

//                 Add your code to emit
//                   - object initializer
//                   - the class methods
//                   - etc...

  if (cgen_debug) cout << "coding class_initObj" << endl;
  code_class_initObj();

  if (cgen_debug) cout << "coding class_method" << endl;
  code_class_method();

}


CgenNodeP CgenClassTable::root()
{
   return probe(Object);
}


///////////////////////////////////////////////////////////////////////
//
// CgenNode methods
//
///////////////////////////////////////////////////////////////////////

CgenNode::CgenNode(Class_ nd, Basicness bstatus, CgenClassTableP ct) :
   class__class((const class__class &) *nd),
   parentnd(NULL),
   children(NULL),
   basic_status(bstatus)
{ 
   stringtable.add_string(name->get_string());          // Add class name to string table
}

void CgenNode::code_classname(ostream &s)
{
	s << WORD; stringtable.lookup_string(name->get_string())->code_ref(s); s << endl;
}

void CgenNode::code_classobj(ostream &s)
{
	s << WORD; emit_protobj_ref(name, s); s << endl;
	s << WORD; emit_init_ref(name, s); s << endl;
}

void CgenNode::code_dispTab(ostream &s)
{ 
	emit_disptable_ref(name, s); s << LABEL;

	std::map<Symbol, Symbol> method2class;
	std::vector<Symbol> methods;
	get_features(Feature_class::METHOD, this, methods, method2class);
	for (int i = 0; i < static_cast<int>(methods.size()); ++i) {
		s << WORD; emit_method_ref(method2class[methods[i]], methods[i], s); s << endl;
	}
}

void CgenNode::code_protObj(ostream &s, int idx)
{
 	s << WORD << "-1" << endl;

	emit_protobj_ref(name, s); s << LABEL;
	s << WORD << idx << endl;

	std::map<Symbol, Symbol> attr2type;
	std::vector<Symbol> attrs;

	get_features(Feature_class::ATTR, this, attrs, attr2type);
	int len = 3 + attrs.size(); 

	s << WORD << len << endl;

	s << WORD; emit_disptable_ref(name, s); s << endl;

	for (int i = 0; i < static_cast<int>(attrs.size()); i++) {
		Symbol type = attr2type[attrs[i]];
		if (Int == type) {
			s << WORD; inttable.lookup_string("0")->code_ref(s); s << endl;
		} else if (Str == type) {
			s << WORD; stringtable.lookup_string("")->code_ref(s); s << endl;
		} else if (Bool == type) {
			s << WORD; falsebool.code_ref(s); s << endl;
		} else {
			s << WORD << 0 << endl;
		}
	}
}


void CgenNode::get_features(int type, CgenNodeP node, std::vector<Symbol> &names, std::map<Symbol, Symbol> &name2class)
{
	if (node->name == No_class)
		return;
	get_features(type, node->parentnd, names, name2class);
	Features features = node->features;
	for (int i = features->first(); features->more(i); i = features->next(i)) {
		Feature feature = features->nth(i);
		if (feature->get_type() == type) {
			Symbol type_decl;
			Symbol name = feature->get_name();
			if (name2class.find(name) == name2class.end()) {
				feature->set_offset(names.size(), type_decl);
				names.push_back(name);
			} else {
				for (int j = 0; j < static_cast<int>(names.size()); j++) {
					if (names[j] == name) {
						feature->set_offset(j, type_decl);
						break;
					}
				}
			}
			if (type == Feature_class::METHOD) {
				name2class[name] = node->name;
			} else if (type == Feature_class::ATTR) {
				name2class[name] = type_decl;
			}
		}
	}
}

void attr_class::set_offset(int n, Symbol &t)
{
	t = type_decl;
	offset = n;
}

void method_class::set_offset(int n, Symbol &_t)
{
	offset = n;
}

int attr_class::get_type()
{
	return Feature_class::ATTR;
}

int method_class::get_type()
{
	return Feature_class::METHOD;
}

void enterfunc(ostream &s)
{
	emit_addiu(SP, SP, -12, s); 
	emit_store(FP, 3, SP, s);
	emit_store(SELF, 2, SP, s);
	emit_store(RA, 1, SP, s);

	emit_addiu(FP, SP, 4, s);

	emit_move(SELF, ACC, s);

	sp_level = 0;

}

void exitfunc(ostream &s, int n)
{
	emit_load(FP, 3, SP, s);
	emit_load(SELF, 2, SP, s);
	emit_load(RA, 1, SP, s);
	emit_addiu(SP, SP, 12 + n * 4, s);
	emit_return(s);
}

void CgenNode::code_initObj(ostream &s)
{
	emit_init_ref(name, s); s << LABEL;
	enterfunc(s);

	if (parentnd->name != No_class) {
		std::ostringstream strout;
		emit_init_ref(parentnd->name, strout);
		emit_jal(strout.str().c_str(), s);
	}

	self_ = this;

	for (int i = features->first(); features->more(i); i = features->next(i)) {
		features->nth(i)->initialization(s);
	}

	emit_move(ACC, SELF, s);

	exitfunc(s, 0);
}



void attr_class::initialization(ostream &s)
{
	if (!init->is_no_expr()) {
		init->code(s);
		emit_store(ACC, 3 + offset, SELF, s);
	}
}

void method_class::emit_method(ostream &s)
{
	emit_method_ref(self_->name, name, s); s << LABEL;

	locations->enterscope();

	int idx = formals->len();
	for (int i = formals->first(); formals->more(i); i = formals->next(i)) {
		locations->addid(formals->nth(i)->get_name(), new int(3 + --idx));
	}

	enterfunc(s);

	expr->code(s);

	exitfunc(s, formals->len());

	locations->exitscope();
}

void CgenNode::code_method(ostream &s)
{
	if (!basic()) {
		self_ = this;

		for (int i = features->first(); features->more(i); i = features->next(i)) {
			features->nth(i)->emit_method(s);
		}
	}
}

int CgenNode::get_offset(Symbol name, int type)
{
	CgenNodeP p = this;
	while (p->name != No_class) {
		Features f = p->features;
		for (int i = f->first(); f->more(i); i = f->next(i)) {
			Feature f0 = f->nth(i);
			if (f0->get_type() == type && f0->get_name() == name) {
				return f0->get_offset();
			}
		}
		p = p->parentnd;
	}
	// never to this
	assert(false);
	return -1;
}

Symbol CgenNode::get_low_method_type(Symbol name)
{

	CgenNodeP p = this;
	while (p->name != No_class) {
		Features f = p->features;
		for (int i = f->first(); f->more(i); i = f->next(i)) {
			Feature f0 = f->nth(i);
			if (f0->get_type() == Feature_class::METHOD && f0->get_name() == name) {
				return p->name;
			}
		}
		p = p->parentnd;
	}
	// never to this
	assert(false);
	return NULL;

}

std::pair<const char *, int> get_location(Symbol name)
{
	int *p = locations->lookup(name);
	if (p != NULL) {
		return std::make_pair(FP, *p);
	} else {
		return std::make_pair(SELF, 3 + self_->get_offset(name, Feature_class::ATTR));
	}
}


/*
void method_class::code(ostream &s, CgenNodeP node)
{
	emit_method_ref(node->name, name, s); s << LABEL;

}

void attr_class::code(ostream &s)
{
	// keep empty
}
*/


//******************************************************************
//
//   Fill in the following methods to produce code for the
//   appropriate expression.  You may add or remove parameters
//   as you wish, but if you do, remember to change the parameters
//   of the declarations in `cool-tree.h'  Sample code for
//   constant integers, strings, and booleans are provided.
//
//*****************************************************************

void call_method(Symbol type_name, Symbol name, ostream &s)
{
	std::ostringstream strout;
	emit_method_ref(type_name, name, strout);
	emit_jal(strout.str().c_str(), s);
}

void load_protobj(Symbol type, ostream &s)
{
	std::ostringstream strout;
	emit_protobj_ref(type, strout);
	emit_load_address(ACC, strout.str().c_str(), s);
}


int get_label()
{
	return LABEL_COUNT++;
}

void assign_class::code(ostream &s) 
{
	expr->code(s);
	std::pair<const char *, int> r = get_location(name);
	emit_store(ACC, r.second, r.first, s);
	/*
	emit_addiu(A1, r.first, r.second * 4, s);
	emit_gc_assign(s);
	*/
}

void emit_void_abort(int line_number, const char *func, ostream &s)
{
	int label = get_label();
	emit_bne(ACC, ZERO, label, s);

	emit_load_string(ACC, stringtable.lookup_string(self_->filename->get_string()), s);
	emit_load_imm(T1, line_number, s);
	emit_jal(func, s);

	emit_label_def(label, s);
}

void code_dispatch(bool is_dynamic, Expression expr, 
		Symbol type_name, Symbol name, Expressions actual, ostream &s)
{
	if (type_name == SELF_TYPE) {
		type_name = self_->name;
	}

	for (int i = actual->first(); actual->more(i); i = actual->next(i)) {
		actual->nth(i)->code(s);
		emit_push(ACC, s);
		sp_level++;
	}

	expr->code(s);
	
	//cout << name << " " << expr->get_type() << endl;

	emit_void_abort(expr->get_line_number(), "_dispatch_abort", s);


	if (is_dynamic) {
		emit_load(T1, 2, ACC, s);
		int offset = class_table->probe(type_name)->get_offset(name, Feature_class::METHOD);
		//cout << name << " " << offset << endl;
		emit_load(T1, offset, T1, s);
		emit_jalr(T1, s);
	} else {
		call_method(class_table->probe(type_name)->get_low_method_type(name), name, s);
	}
	sp_level -= actual->len();
}

void static_dispatch_class::code(ostream &s) 
{
	code_dispatch(false, expr, type_name, name, actual, s);
}

void dispatch_class::code(ostream &s) 
{
	code_dispatch(true, expr, expr->get_type(), name, actual, s);
}


void cond_class::code(ostream &s) 
{
	pred->code(s);
	int l1 = get_label(), l2 = get_label();
	emit_load(ACC, 3, ACC, s);
	emit_beqz(ACC, l1, s);
	then_exp->code(s);
	emit_branch(l2, s);
	emit_label_def(l1, s);
	else_exp->code(s);
	emit_label_def(l2, s);
}

void loop_class::code(ostream &s) 
{
	int loop_label  = get_label(), outloop_label = get_label();
	emit_label_def(loop_label, s);
	pred->code(s);
	emit_load(ACC, 3, ACC, s);
	emit_beqz(ACC, outloop_label, s);

	body->code(s);

	emit_branch(loop_label, s);
	emit_label_def(outloop_label, s);

}

void beq_tag(Symbol x, std::vector<int> &labels, ostream &s)
{
	emit_load_imm(T2, n2tag[x], s);
	int label = get_label();
	labels.push_back(label);
	emit_beq(T1, T2, label, s);
}

void typcase_class::code(ostream &s) 
{
	expr->code(s);

	emit_void_abort(expr->get_line_number(), "_case_abort2", s);

	// load tag
	emit_load(T1, 0, ACC, s);

	int loop_label = get_label();

	emit_label_def(loop_label, s);

	std::vector<int> labels;
	for (int i = cases->first(); cases->more(i); i = cases->next(i)) {
		branch_class *b = static_cast<branch_class *>(cases->nth(i));
		beq_tag(b->type_decl, labels, s);
	}
	emit_load_address(T2, CLASSINHERITTAB, s);
	emit_sll(T1, T1, 2, s);
	emit_addu(T2, T2, T1, s);
	emit_load(T1, 0, T2, s);
	emit_load_imm(T2, n2tag[No_class], s);
	emit_bne(T1, T2, loop_label, s);

	emit_move(ACC, SELF, s);
	emit_jal("_case_abort", s);


	int final_label = get_label();

	int idx = 0;
	for (int i = cases->first(); cases->more(i); i = cases->next(i)) {
		branch_class *b = static_cast<branch_class *>(cases->nth(i));
		emit_label_def(labels[idx++], s);

		emit_push(ACC, s);
		sp_level++;

		locations->enterscope();
		locations->addid(b->name, new int(-sp_level));

		b->expr->code(s);

		locations->exitscope();

		emit_addiu(SP, SP, 4, s);
		sp_level--;

		emit_branch(final_label, s);
	}

	emit_label_def(final_label, s);


}

void block_class::code(ostream &s) 
{
	for (int i = body->first(); body->more(i); i = body->next(i)) {
		body->nth(i)->code(s);
	}
}

void let_class::code(ostream &s) 
{
	init->code(s);


	if (init->is_no_expr()) {
		std::ostringstream strout;
		if (Int == type_decl) {
			inttable.lookup_string("0")->code_ref(strout);
		} else if (Str == type_decl) {
			stringtable.lookup_string("")->code_ref(strout); 
		} else if (Bool == type_decl) {
			falsebool.code_ref(strout); 
		} 			

		if (strout.str().length() != 0) {
			emit_load_address(ACC, strout.str().c_str(), s);
		} else {
			emit_load_imm(ACC, 0, s);
		}
	} 
	emit_push(ACC, s);
	sp_level++;

	locations->enterscope();

	locations->addid(identifier, new int(-sp_level));

	body->code(s);

	locations->exitscope();

	emit_addiu(SP, SP, 4, s);

	sp_level--;
}

void save_t1(ostream &s)
{
	emit_push(T1, s);
	sp_level++;
}

void load_t1(ostream &s)
{
	emit_addiu(SP, SP, 4, s);
	sp_level--;
	emit_load(T1, 0, SP, s);
}

void code_2expr(Expression e1, Expression e2, ostream &s)
{
	e1->code(s);
	emit_fetch_int(T1, ACC, s);
	save_t1(s);
	e2->code(s);
	load_t1(s);
	emit_fetch_int(ACC, ACC, s);
}


void code_arith(Expression e1, Expression e2, void (op)(char *, char *, char *, ostream&), ostream &s)
{
	code_2expr(e1, e2, s);
	op(T1, T1, ACC, s);

	load_protobj(Int, s);
	save_t1(s);
	call_method(Object, ::copy, s);
	load_t1(s);

	emit_store_int(T1, ACC, s);
}

void code_compare(Expression e1, Expression e2, void (op)(char *, char *, int, ostream &), ostream &s)
{
	code_2expr(e1, e2, s);
	emit_load_bool(T2, truebool, s);
	int label = get_label();
	op(T1, ACC, label, s);
	emit_load_bool(T2, falsebool, s);
	emit_label_def(label, s);
	emit_move(ACC, T2, s);
}


void plus_class::code(ostream &s) 
{
	code_arith(e1, e2, emit_add, s);
}

void sub_class::code(ostream &s) 
{
	code_arith(e1, e2, emit_sub, s);
}

void mul_class::code(ostream &s) 
{
	code_arith(e1, e2, emit_mul, s);
}

void divide_class::code(ostream &s) 
{
	code_arith(e1, e2, emit_div, s);
}

void neg_class::code(ostream &s) 
{
	e1->code(s);
	call_method(Object, ::copy, s);
	emit_fetch_int(T1, ACC, s);
	emit_neg(T1, T1, s);
	emit_store_int(T1, ACC, s);
}

void lt_class::code(ostream &s) 
{
	code_compare(e1, e2, emit_blt, s);
}

void eq_class::code(ostream &s) 
{
	e1->code(s);
	emit_move(T1, ACC, s);
	save_t1(s);
	e2->code(s);
	load_t1(s);
	emit_move(T2, ACC, s);

	emit_load_bool(ACC, truebool, s);

	int label = get_label();

	emit_beq(T1, T2, label, s);

	emit_load_bool(A1, falsebool, s);
	emit_jal("equality_test", s);

	emit_label_def(label, s);
	
}

void leq_class::code(ostream &s) 
{
	code_compare(e1, e2, emit_bleq, s);
}

void comp_class::code(ostream &s) 
{
	e1->code(s);
	emit_load(T1, DEFAULT_OBJFIELDS, ACC, s);
	emit_load_bool(ACC, truebool, s);
	int label = get_label();
	emit_beqz(T1, label, s);
	emit_load_bool(ACC, falsebool, s);
	emit_label_def(label, s);
}

void int_const_class::code(ostream& s)  
{
  //
  // Need to be sure we have an IntEntry *, not an arbitrary Symbol
  //
  emit_load_int(ACC,inttable.lookup_string(token->get_string()),s);
}

void string_const_class::code(ostream& s)
{
  emit_load_string(ACC,stringtable.lookup_string(token->get_string()),s);
}

void bool_const_class::code(ostream& s)
{
  emit_load_bool(ACC, BoolConst(val), s);
}

void code_obj_offset(int id, ostream& s)
{
	emit_load(T1, 0, SELF, s);
	emit_sll(T1, T1, 3, s);
	emit_addiu(T1, T1, id * 4, s);
	emit_load_address(T2, CLASSOBJTAB, s);
	emit_addu(T2, T2, T1, s);
}

void new__class::code(ostream &s) 
{
	Symbol type = type_name;
	if (type == SELF_TYPE) {
		code_obj_offset(0, s);
		emit_load(ACC, 0, T2, s);
		call_method(Object, ::copy, s);
		code_obj_offset(1, s);
		emit_load(T2, 0, T2, s);
		emit_jalr(T2, s);
	} else {
		load_protobj(type, s);
		call_method(Object, ::copy, s);
		std::ostringstream strout;
		emit_init_ref(type, strout);
		emit_jal(strout.str().c_str(), s);
	}
}

void isvoid_class::code(ostream &s) 
{
	e1->code(s);
	emit_move(T1, ACC, s);
	emit_load_bool(ACC, truebool, s);
	int label = get_label(); 
	emit_beqz(T1, label, s);
	emit_load_bool(ACC, falsebool, s);
	emit_label_def(label, s);
}

void no_expr_class::code(ostream &s) 
{
	// keep empty
}

void object_class::code(ostream &s) 
{
	if (name == self) {
		emit_move(ACC, SELF, s);
	} else {
		std::pair<const char *, int> ret = get_location(name);
		emit_load(ACC, ret.second, ret.first, s);
	}
}




