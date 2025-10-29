// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME EventReaderDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "EventReader.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_EventReader(void *p = nullptr);
   static void *newArray_EventReader(Long_t size, void *p);
   static void delete_EventReader(void *p);
   static void deleteArray_EventReader(void *p);
   static void destruct_EventReader(void *p);
   static void streamer_EventReader(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventReader*)
   {
      ::EventReader *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EventReader >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("EventReader", ::EventReader::Class_Version(), "EventReader.h", 18,
                  typeid(::EventReader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EventReader::Dictionary, isa_proxy, 16,
                  sizeof(::EventReader) );
      instance.SetNew(&new_EventReader);
      instance.SetNewArray(&newArray_EventReader);
      instance.SetDelete(&delete_EventReader);
      instance.SetDeleteArray(&deleteArray_EventReader);
      instance.SetDestructor(&destruct_EventReader);
      instance.SetStreamerFunc(&streamer_EventReader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EventReader*)
   {
      return GenerateInitInstanceLocal(static_cast<::EventReader*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::EventReader*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr EventReader::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *EventReader::Class_Name()
{
   return "EventReader";
}

//______________________________________________________________________________
const char *EventReader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventReader*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int EventReader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventReader*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EventReader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventReader*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EventReader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventReader*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void EventReader::Streamer(TBuffer &R__b)
{
   // Stream an object of class EventReader.

   TSelector::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EventReader(void *p) {
      return  p ? new(p) ::EventReader : new ::EventReader;
   }
   static void *newArray_EventReader(Long_t nElements, void *p) {
      return p ? new(p) ::EventReader[nElements] : new ::EventReader[nElements];
   }
   // Wrapper around operator delete
   static void delete_EventReader(void *p) {
      delete (static_cast<::EventReader*>(p));
   }
   static void deleteArray_EventReader(void *p) {
      delete [] (static_cast<::EventReader*>(p));
   }
   static void destruct_EventReader(void *p) {
      typedef ::EventReader current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_EventReader(TBuffer &buf, void *obj) {
      ((::EventReader*)obj)->::EventReader::Streamer(buf);
   }
} // end of namespace ROOT for class ::EventReader

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_EventReaderDict() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_EventReaderDict) = RecordReadRules_EventReaderDict();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_EventReaderDict));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_EventReaderDict_Impl() {
    static const char* headers[] = {
"EventReader.h",
nullptr
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/dybfs2/users/lijj16/nHBackEndAnalysis/LS_U_Th_Level/AnalysisLSAlphas/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "EventReaderDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$EventReader.h")))  EventReader;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "EventReaderDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "EventReader.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"EventReader", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("EventReaderDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_EventReaderDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_EventReaderDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_EventReaderDict() {
  TriggerDictionaryInitialization_EventReaderDict_Impl();
}
