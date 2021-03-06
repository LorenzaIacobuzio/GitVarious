#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include<iostream> 

void CopyDir(TDirectory *source) {
  //copy all objects and subdirs of directory source as a subdir of the current directory   
  //source->ls();
  TDirectory *savdir = gDirectory;
  TDirectory *adir = savdir->mkdir(source->GetName());
  adir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir);
      adir->cd();
    }
    else {
      source->cd();
      TObject *obj = key->ReadObj();
      adir->cd();
      obj->Write();
      delete obj;
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}


void CopyDir(TDirectory *source, const char *directory, const char *tree) {
  //copy all objects and subdirs of directory source as a subdir of the current directory   
  TDirectory *savdir = gDirectory;
  TDirectory *adir = gDirectory;
  adir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if(!strcmp(key->GetName(),directory)){
      if (cl->InheritsFrom(TDirectory::Class())) {
	source->cd(key->GetName());
	TDirectory *subdir = gDirectory;
	adir->cd();
	CopyDir(subdir);
	adir->cd();
      }
    }
    else if(!strcmp(key->GetName(),tree)){
      if (cl->InheritsFrom(TTree::Class())) {
	TTree *T = (TTree*)source->Get(key->GetName());
	adir->cd();
	TTree *newT = T->CloneTree(-1,"fast");
	newT->Write();
      }
    }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void CopyDir(TDirectory *source, const char *directory) {
  //copy all objects and subdirs of directory source as a subdir of the current directory   
  TDirectory *savdir = gDirectory;
  TDirectory *adir = gDirectory;
  adir->cd();
  //loop on all entries of this directory
  TKey *key;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
    const char *classname = key->GetClassName();
    TClass *cl = gROOT->GetClass(classname);
    if (!cl) continue;
    if(strcmp(key->GetName(),directory)) continue;
    if (cl->InheritsFrom(TDirectory::Class())) {
      source->cd(key->GetName());
      TDirectory *subdir = gDirectory;
      adir->cd();
      CopyDir(subdir);
      adir->cd();
    } 
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}

void CopyFile(const char *fname, const char *dir, const char *tree) {
   //Copy all objects and subdirs of file fname as a subdir of the current directory
   TDirectory *target = gDirectory;
   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) {
      printf("Cannot copy file: %s\n",fname);
      target->cd();
      return;
   }
   target->cd();
   CopyDir(f,dir,tree);
   delete f;
   target->cd();
}
void CopyFile(const char *fname, const char *dir) {
   //Copy all objects and subdirs of file fname as a subdir of the current directory
   TDirectory *target = gDirectory;
   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) {
      printf("Cannot copy file: %s\n",fname);
      target->cd();
      return;
   }
   target->cd();
   CopyDir(f,dir);
   delete f;
   target->cd();
}
int main (int argc, char **argv) {
  if(argc<4){
    std::cout<<"Wrong number of parameters: ./copyFile input.root output.root directory1 (directory2)"<<std::endl;
    exit(EXIT_FAILURE);
  } 
  TFile *f = new TFile(argv[2],"recreate");
  if(argc==4){
    CopyFile(argv[1],argv[3]);
  }
  if(argc==5){
    CopyFile(argv[1],argv[3],argv[4]);
  }
  delete f;
  
  return 0;
}
