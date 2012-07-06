/* 
   


 */


extern "C" {
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>
}

/**
Defines:
*/

#define LOCUSLEN         5
#define ALLELELEN        4
#define NONGENOTYPECOLS  6

///names of parameters in the simulation.  Used as names in lists

#define INTEGERPARAMS    "intparam"
#define SWITCHPARAMS     "switchparam"
#define FLOATPARAMS      "floatparam"
#define DEMOPARAMS       "demography"
#define LOCIPARAMS       "loci"
#define INDPARAMS        "individuals"

#define HABNAMES         "habitats"     
#define STAGENAME 	 "stages"       
#define	LNUMNAME	 "locusnum"     
#define	ENUMNAME	 "numepochs"    
#define	CGNAME		 "currentgen"   
#define	CENAME		 "currentepoch" 
#define	FINALAGE	 "totalgens"    
#define	DNUMNAME	 "numdemos"     
#define	MAXLANDNAME	 "maxlandsize"  
#define NEXTIDNAME       "nextid"

#define TYPENAME	 "type"    
#define	PLOIDYNAME	 "ploidy"  
#define	TRANSNAME	 "trans"   
#define	RATENAME	 "rate"    
#define	ALISTNAME 	 "alleles" 

#define AINDXNAME        "aindex"  
#define	ABIRTHNAME	 "birth"   
#define	PROPNAME	 "prop"    
#define	STATENAME	 "state"   

#define LOCALDEMNM       "localdem"
#define LOCALDEMKNM      "localdemK"
#define EPOCHDEMNM       "epochs"

#define LCLSMATNM        "LocalS"
#define LCLRMATNM        "LocalR"
#define LCLMMATNM        "LocalM"

#define RNDCHSNAME       "RndChooseProb"
#define	SGENAME	         "StartGen"     
#define	EXTINCTNAME	 "Extinct"      
#define	CARRYNAME	 "Carry"        
#define	LPNAME		 "Localprob"    
#define SNAME            "S"    
#define	RNAME		 "R"            
#define	MNAME      	 "M"            

#define RANDEPOCHN       "randepoch"
#define RANDDEMON 	 "randdemo" 
#define MULTPNAME        "multp"
///KKM 5.20.05..............................................................
#define DENSDEP     "densdepdemo"
///.........................................................................
#define SELFRATENAME     "selfing"

///input/output of landscapes to/from metasim lib
extern "C" SEXP read_landscape(SEXP fn);
//extern "C" SEXP convert_metasim_to_R(Landscape_statistics &L);
extern "C" SEXP write_landscape(SEXP fn, SEXP Rland);
//extern "C" void convert_R_to_metasim(SEXP Rland, Landscape_statistics &L);
extern "C" SEXP getListElement(SEXP list, const char *str);

///simulations
///run metasim on the landscape a certain number of times
extern "C" SEXP iterate_landscape(SEXP numit, SEXP Rland, SEXP cmpress, SEXP bypop);
///perform survival step on the landscape
extern "C" SEXP survive_landscape(SEXP Rland);
///perform reproduce step on the landscape
extern "C" SEXP reproduce_landscape(SEXP Rland);
///perform carry step on the landscape
extern "C" SEXP carry_landscape(SEXP Rland);
///perform extinct step on the landscape
extern "C" SEXP extinct_landscape(SEXP Rland);
///advance the landscape through time
extern "C" SEXP advance_landscape(SEXP Rland);


extern "C" SEXP populate_Rland(SEXP Rland, SEXP Population_sizes);

extern "C" SEXP clean_landscape(SEXP Rland);
extern "C" SEXP compress_landscape(SEXP Rland);

extern "C" SEXP num_demo_cols();
///utility functions
///convert a landscape into a format that the weir fst calculations in R can use.
extern "C" SEXP l2w(SEXP Rland, SEXP numind);

///
///this is C code to calculate relatedness
///

extern "C" SEXP relateinternal(SEXP ind, SEXP acnp);
 

///output functions
///convert landscapes to various formats
///All take a filename, a landscape, and a number of indiviudals to sample per habitat
extern "C" SEXP writeGDA(SEXP fn, SEXP Rland, SEXP ni);
extern "C" SEXP writeArlequinDip(SEXP fn, SEXP Rland, SEXP ni);
extern "C" SEXP writeArlequinHap(SEXP fn, SEXP Rland, SEXP ni);
extern "C" SEXP writeBIOSYS(SEXP fn, SEXP Rland, SEXP ni);
extern "C" SEXP writeGenPop(SEXP fn, SEXP Rland, SEXP ni);
extern "C" SEXP writeR(SEXP fn, SEXP Rland, SEXP ni);
extern "C" SEXP writeMigrateDip(SEXP fn, SEXP Rland, SEXP ni);
extern "C" SEXP writeReRat(SEXP fn, SEXP Rland, SEXP ni);

extern "C" SEXP test(SEXP mat1, SEXP mat2);

