/*
Copyright (C) [2023] [Alexander Volya]
This code is licensed under the GNU GPL v3, see LICENSE file 
For more details and acknowledgments, see the NOTICE file.
*/


/*! 
\author Alexander Volya  <http://www.volya.net>
\file label.cxx
\brief Control Human-Computer labeling 
*/
/*
Purpose: Organize labeling of file names
2/19/2010 parity needs some fixing in label, PhotonLabel is ok
1/13/06 created
03/13/2006 temporary_name global variable is created
          Function AddExtension is created to be used with file extension
	  New function that read quantum numbers
07/25/06  AddChar(int,char) and (char,int)
I also made PhotonLabel but this is not easy because M is there M1_M2T0N0+
 */
#ifndef _LABEL_CXX_
#define _LABEL_CXX_
#include <cstdio>
using std::strcpy;
using std::strcat;
using std::snprintf;
namespace SM{
  const int default_label_length=256;
  char generic_temporary_name[default_label_length]; //use in SM independent function
  char generic_temporary_name1[default_label_length];
  char temporary_name[default_label_length]; //global char
  char temporary_name1[default_label_length]; //global char
  char lname[]={'s','p','d','f','g','h','i','j','k','l'}; ///<Array of orbital names
  char charmap[]={'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','y','Z'}; ///<ASCII character map
/*! \brief Spherical Single-Particle label
The form of output p0s1 for proton on 0s j=1/2 orbit 
*/
  int SPLabel(char *lbl, ///<Output label
 int tz, ///<isospin projection (proton or neutron)
int main_qn, ///< main quantum number 
int ll, ///<orbital angular momentum
int jj ///<total angular momentum
) {
    //labl must exist already
    char  stmp[10]; //temporary label
    if (tz==1) strcpy (lbl,"n"); 
    else strcpy(lbl, "p");  //first symbol for proton or neutron
    snprintf (stmp, 10, "%u", main_qn);
    strcat (lbl,stmp);
    snprintf (stmp, 10, "%c", lname[ll]); //put simbol for orbital
    strcat (lbl,stmp);
     snprintf (stmp,10, "%u", jj);
    strcat (lbl,stmp);
    return 1;
  }
  /*! \brief Spherical Single-Particle label, isospin invariant
The form of output p0s1 for proton on 0s j=1/2 orbit 
*/
  int SPLabel(char *lbl, ///<Output label
int main_qn, ///< main quantum number 
int ll, ///<orbital angular momentum
int jj ///<total angular momentum
) {
    //labl must exist already
    char  stmp[10]; //temporary label
    strcpy (lbl,""); 
	snprintf (stmp,10, "%u", main_qn);
    strcat (lbl,stmp);
    snprintf (stmp, 10,"%c", lname[ll]); //put simbol for orbital
    strcat (lbl,stmp);
     snprintf (stmp,10, "%u", jj);
    strcat (lbl,stmp);
    return 1;
  }
  
/*!\brief Spherical Single-Particle label, return pointer
*/
 char* SPLabel(int tz, int main_qn, int ll, int jj) {
   SPLabel(temporary_name, tz, main_qn, ll, jj);
   return temporary_name;
  }
  /*!\brief Spherical Single-Particle label, return pointer, isospininvariant
*/
 char* SPLabel(int main_qn, int ll, int jj) {
   SPLabel(temporary_name, main_qn, ll, jj);
   return temporary_name;
  }
/*!\brief Add Spherical Single-Particle label to other label
*/
  int AppendSPLabel(char *lbl, int tz, int main_qn, int ll, int jj) {
    //labl must exist already
    char  stmp[20]; //temporary label
    SPLabel(stmp, tz,  main_qn, ll, jj);
    strcat (lbl, "_");
    strcat (lbl,stmp);
    return 1;
  }
/*!\brief Read Spherical Single-Particle label
*/
  int ReadSPLabel(char *lbl, int &tz, int &main_qn, int &ll, int &jj){
   const int czero=int('0');
    if (lbl[0]=='n') {tz=1;}
    else {tz=-1;}
    main_qn=0;
    int i=0;
     for (;lbl[i]!='\0';i++) 
       	{ 
	  i++; //enetring value
	  // if (lbl[i]=='-') {signflag=true; i++;}
	  while ((lbl[i]>=czero)&&(lbl[i]<(czero+10))) 
	    {main_qn*=10; main_qn+=lbl[i]-czero; i++; } //while number read it
	  break; 
	}
	for (ll=0;ll<8;ll++) if (lbl[i]==lname[ll]) break;
	jj=0;
      for (;lbl[i]!='\0';i++) {
       	{ 
	  i++; //enetring value
	  // if (lbl[i]=='-') {signflag=true; i++;}
	  while ((lbl[i]>=czero)&&(lbl[i]<(czero+10))) 
	    {jj*=10; jj+=lbl[i]-czero; i++; } //while number read it
	  break; 
	}
    
      }
return 0;
  }
/*!\brief add two char* strings
*/
  char * AddExtension(const char *label, const char *suffix){
    strcpy (generic_temporary_name, label); //put label
    strcat (generic_temporary_name, suffix);
    return generic_temporary_name;
  }
/*!\brief add three char* strings
*/
  char * AddChar(const char *a, const char *b, const char *c){
    strcpy (generic_temporary_name, a); //put label
    strcat (generic_temporary_name, b);
    strcat (generic_temporary_name, c);
    return generic_temporary_name;
  }
/*!\brief add two char* strings
*/
  char * AddChar(const char *a, const char *b){
    strcpy (generic_temporary_name, a); //put label
    strcat (generic_temporary_name, b);
    return generic_temporary_name;
  }
/*!\brief add integer to char* at end
*/
    char * AddChar(const char *a, const int b){
        strcpy (generic_temporary_name, a); //put label
      snprintf (generic_temporary_name1,default_label_length, "%i",b); //put b into temporary name
      strcat (generic_temporary_name, generic_temporary_name1);
    return generic_temporary_name;
  }
/*!\brief add integer to char* at begin
*/
     char * AddChar(const int b, const char* a){
	  static char temporary_char[512];
      snprintf (temporary_char,512, "%i",b); //put b into temporary name
      strcat (temporary_char, a);
    return temporary_char;
  }

#ifdef __QN__CXX__
  int QNLabel(char *lbl, QN &q) {
    if(q.P==-1) snprintf (lbl,default_label_length, "M%iT%iN%i-", int(q.Jz),int(q.Tz),int(q.N));
    else snprintf (lbl,default_label_length, "M%iT%iN%i+", int(q.Jz),int(q.Tz),int(q.N));
      //sprintf (lbl, "%c", char(q.Jz));
      //sprintf (lbl, "%c", char(q.Tz)); 
    return 1;
  }
  char * QNLabel(QN &q) {
    QNLabel(temporary_name, q);
    return temporary_name;
  }

 
   int AppendQNLabel(char *lbl, QN &q) {
     QNLabel(temporary_name, q);
     strcat(lbl, temporary_name);
    return 1;
  }

   char * PhotonLabel(QN &q) {
     //electic or magnetic?
     if (((q.J>>1)&1)^q.P) { //natural parity check
       //magnetic
       strcpy(temporary_name1,AddChar("M",(q.J/2)));
      strcat (temporary_name1, "_");
       AppendQNLabel(temporary_name1, q);
       return temporary_name1;
     }
     else 
	 //electric
       strcpy(temporary_name1,AddChar("E",(q.J/2)));
       strcat (temporary_name1, "_");
       AppendQNLabel(temporary_name1, q);
       return temporary_name1;
 
  }


  int ReadintLabel(char *lbl, char signal) {
    //null-terminated string
     const int czero=int('0');
     int iret=0;
     int signflag=false;
     for (int i=0;lbl[i]!='\0';i++) {
     if (lbl[i]==signal) //read moment
	{ 
	  i++; //enetring value
	    if (lbl[i]=='-') {signflag=true; i++;}
	  while ((lbl[i]>=czero)&&(lbl[i]<(czero+10))) 
	    {iret*=10; iret+=lbl[i]-czero; i++; } //while number read it
	  break; 
	}
     }
     if (signflag) iret=-iret;
     return iret;
  }
  int ReadQNLabel(char *lbl, QN &Q) {
    Q.Jz=ReadintLabel(lbl,'M');
    Q.Tz=ReadintLabel(lbl,'T');
    Q.N=ReadintLabel(lbl,'N');
    for(int i=0;;i++) if(lbl[i]=='\0') { 
    if (lbl[i-1]=='-') Q.P=-1; else Q.P=1;
    break ;
    }
    return 1;
  }

#endif
#ifdef __SPACE__CXX__
  int SpaceLabel(char *lbl, Space &q) {
 if (q.type==1) strcpy (lbl,"f"); else strcpy(lbl,"b");
 //char  stmp[256];
 for (int i=0;i<q.levels;i++) {
 std::snprintf (temporary_name, default_label_length,"%c", q.j[i]);
 strcat (lbl,temporary_name);
 }
 return 1;
  }
  char * SpaceLabel(Space &q) {
    SpaceLabel(temporary_name, q);
    return temporary_name;
  }
  int AppendSpaceLabel(char *lbl, Space &q) {
     SpaceLabel(temporary_name, q);
     strcat(lbl, temporary_name);
    return 1;
  }
#endif
  
}
#endif

