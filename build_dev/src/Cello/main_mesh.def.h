

/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_exit_2_closure : public SDAG::Closure {
            int count_blocks;


      p_exit_2_closure() {
        init();
      }
      p_exit_2_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return count_blocks;}
      void pup(PUP::er& __p) {
        __p | count_blocks;
        packClosure(__p);
      }
      virtual ~p_exit_2_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_exit_2_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_checkpoint_output_3_closure : public SDAG::Closure {
            int count;
            std::string dir;


      p_checkpoint_output_3_closure() {
        init();
      }
      p_checkpoint_output_3_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return count;}
            std::string & getP1() { return dir;}
      void pup(PUP::er& __p) {
        __p | count;
        __p | dir;
        packClosure(__p);
      }
      virtual ~p_checkpoint_output_3_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_checkpoint_output_3_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_initial_exit_4_closure : public SDAG::Closure {
      

      p_initial_exit_4_closure() {
        init();
      }
      p_initial_exit_4_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_initial_exit_4_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_initial_exit_4_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_adapt_enter_5_closure : public SDAG::Closure {
      

      p_adapt_enter_5_closure() {
        init();
      }
      p_adapt_enter_5_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_enter_5_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_enter_5_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_adapt_called_6_closure : public SDAG::Closure {
      

      p_adapt_called_6_closure() {
        init();
      }
      p_adapt_called_6_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_called_6_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_called_6_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_adapt_end_7_closure : public SDAG::Closure {
      

      p_adapt_end_7_closure() {
        init();
      }
      p_adapt_end_7_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_end_7_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_end_7_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_adapt_update_8_closure : public SDAG::Closure {
      

      p_adapt_update_8_closure() {
        init();
      }
      p_adapt_update_8_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_update_8_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_update_8_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_adapt_exit_9_closure : public SDAG::Closure {
      

      p_adapt_exit_9_closure() {
        init();
      }
      p_adapt_exit_9_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_exit_9_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_exit_9_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_compute_enter_10_closure : public SDAG::Closure {
      

      p_compute_enter_10_closure() {
        init();
      }
      p_compute_enter_10_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_compute_enter_10_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_compute_enter_10_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_compute_continue_11_closure : public SDAG::Closure {
      

      p_compute_continue_11_closure() {
        init();
      }
      p_compute_continue_11_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_compute_continue_11_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_compute_continue_11_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_compute_exit_12_closure : public SDAG::Closure {
      

      p_compute_exit_12_closure() {
        init();
      }
      p_compute_exit_12_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_compute_exit_12_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_compute_exit_12_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_output_enter_13_closure : public SDAG::Closure {
      

      p_output_enter_13_closure() {
        init();
      }
      p_output_enter_13_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_output_enter_13_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_output_enter_13_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_output_exit_14_closure : public SDAG::Closure {
      

      p_output_exit_14_closure() {
        init();
      }
      p_output_exit_14_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_output_exit_14_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_output_exit_14_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_stopping_enter_15_closure : public SDAG::Closure {
      

      p_stopping_enter_15_closure() {
        init();
      }
      p_stopping_enter_15_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_stopping_enter_15_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_stopping_enter_15_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_stopping_balance_16_closure : public SDAG::Closure {
      

      p_stopping_balance_16_closure() {
        init();
      }
      p_stopping_balance_16_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_stopping_balance_16_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_stopping_balance_16_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_stopping_exit_17_closure : public SDAG::Closure {
      

      p_stopping_exit_17_closure() {
        init();
      }
      p_stopping_exit_17_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_stopping_exit_17_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_stopping_exit_17_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_text_file_write_18_closure : public SDAG::Closure {
            int nd;
            char *dir;
            int nf;
            char *file;
            int nl;
            char *line;
            int count;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_text_file_write_18_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_text_file_write_18_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return nd;}
            char *& getP1() { return dir;}
            int & getP2() { return nf;}
            char *& getP3() { return file;}
            int & getP4() { return nl;}
            char *& getP5() { return line;}
            int & getP6() { return count;}
      void pup(PUP::er& __p) {
        __p | nd;
        __p | nf;
        __p | nl;
        __p | count;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> nd;
  implP|nd;
  int impl_off_dir, impl_cnt_dir;
  implP|impl_off_dir;
  implP|impl_cnt_dir;
  PUP::detail::TemporaryObjectHolder<int> nf;
  implP|nf;
  int impl_off_file, impl_cnt_file;
  implP|impl_off_file;
  implP|impl_cnt_file;
  PUP::detail::TemporaryObjectHolder<int> nl;
  implP|nl;
  int impl_off_line, impl_cnt_line;
  implP|impl_off_line;
  implP|impl_cnt_line;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
          impl_buf+=CK_ALIGN(implP.size(),16);
          dir = (char *)(impl_buf+impl_off_dir);
          file = (char *)(impl_buf+impl_off_file);
          line = (char *)(impl_buf+impl_off_line);
        }
      }
      virtual ~p_text_file_write_18_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_text_file_write_18_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Main::p_exit_19_closure : public SDAG::Closure {
      

      p_exit_19_closure() {
        init();
      }
      p_exit_19_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_exit_19_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_exit_19_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */



/* DEFS: readonly CProxy_Main proxy_main;
 */
extern CProxy_Main proxy_main;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_proxy_main(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|proxy_main;
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
void p_exit(int count_blocks);
void p_checkpoint_output(int count, const std::string &dir);
void p_initial_exit();
void p_adapt_enter();
void p_adapt_called();
void p_adapt_end();
void p_adapt_update();
void p_adapt_exit();
void p_compute_enter();
void p_compute_continue();
void p_compute_exit();
void p_output_enter();
void p_output_exit();
void p_stopping_enter();
void p_stopping_balance();
void p_stopping_exit();
void p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count);
void p_exit();
Main(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Main::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Main(CkArgMsg* impl_msg);
 */
CkChareID CProxy_Main::ckNew(CkArgMsg* impl_msg, int impl_onPE)
{
  CkChareID impl_ret;
  CkCreateChare(CkIndex_Main::__idx, CkIndex_Main::idx_Main_CkArgMsg(), impl_msg, &impl_ret, impl_onPE);
  return impl_ret;
}
void CProxy_Main::ckNew(CkArgMsg* impl_msg, CkChareID* pcid, int impl_onPE)
{
  CkCreateChare(CkIndex_Main::__idx, CkIndex_Main::idx_Main_CkArgMsg(), impl_msg, pcid, impl_onPE);
}
  CProxy_Main::CProxy_Main(CkArgMsg* impl_msg, int impl_onPE)
{
  CkChareID impl_ret;
  CkCreateChare(CkIndex_Main::__idx, CkIndex_Main::idx_Main_CkArgMsg(), impl_msg, &impl_ret, impl_onPE);
  ckSetChareID(impl_ret);
}

// Entry point registration function
int CkIndex_Main::reg_Main_CkArgMsg() {
  int epidx = CkRegisterEp("Main(CkArgMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_Main_CkArgMsg), CMessage_CkArgMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkArgMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Main::_call_Main_CkArgMsg(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  new (impl_obj_void) Main((CkArgMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_exit(int count_blocks);
 */
void CProxy_Main::p_exit(int count_blocks, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int count_blocks
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count_blocks;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count_blocks;
  }
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_exit_marshall2(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_exit_marshall2(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_exit_marshall2(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_exit_marshall2() {
  int epidx = CkRegisterEp("p_exit(int count_blocks)",
      reinterpret_cast<CkCallFnPtr>(_call_p_exit_marshall2), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_exit_marshall2);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_exit_marshall2);

  return epidx;
}

void CkIndex_Main::_call_p_exit_marshall2(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int count_blocks*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count_blocks;
  implP|count_blocks;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_exit(std::move(count_blocks.t));
}
int CkIndex_Main::_callmarshall_p_exit_marshall2(char* impl_buf, void* impl_obj_void) {
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int count_blocks*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count_blocks;
  implP|count_blocks;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_exit(std::move(count_blocks.t));
  return implP.size();
}
void CkIndex_Main::_marshallmessagepup_p_exit_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int count_blocks*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count_blocks;
  implP|count_blocks;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("count_blocks");
  implDestP|count_blocks;
}
PUPable_def(SINGLE_ARG(Closure_Main::p_exit_2_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_checkpoint_output(int count, const std::string &dir);
 */
void CProxy_Main::p_checkpoint_output(int count, const std::string &dir, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int count, const std::string &dir
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
  }
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_checkpoint_output_marshall3(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_checkpoint_output_marshall3(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_checkpoint_output_marshall3(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_checkpoint_output_marshall3() {
  int epidx = CkRegisterEp("p_checkpoint_output(int count, const std::string &dir)",
      reinterpret_cast<CkCallFnPtr>(_call_p_checkpoint_output_marshall3), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_checkpoint_output_marshall3);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_checkpoint_output_marshall3);

  return epidx;
}

void CkIndex_Main::_call_p_checkpoint_output_marshall3(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int count, const std::string &dir*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  PUP::detail::TemporaryObjectHolder<std::string> dir;
  implP|dir;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_checkpoint_output(std::move(count.t), std::move(dir.t));
}
int CkIndex_Main::_callmarshall_p_checkpoint_output_marshall3(char* impl_buf, void* impl_obj_void) {
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int count, const std::string &dir*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  PUP::detail::TemporaryObjectHolder<std::string> dir;
  implP|dir;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_checkpoint_output(std::move(count.t), std::move(dir.t));
  return implP.size();
}
void CkIndex_Main::_marshallmessagepup_p_checkpoint_output_marshall3(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int count, const std::string &dir*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  PUP::detail::TemporaryObjectHolder<std::string> dir;
  implP|dir;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("count");
  implDestP|count;
  if (implDestP.hasComments()) implDestP.comment("dir");
  implDestP|dir;
}
PUPable_def(SINGLE_ARG(Closure_Main::p_checkpoint_output_3_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_initial_exit();
 */
void CProxy_Main::p_initial_exit(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_initial_exit_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_initial_exit_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_initial_exit_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_initial_exit_void() {
  int epidx = CkRegisterEp("p_initial_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_initial_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_initial_exit_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_initial_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_initial_exit_4_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_enter();
 */
void CProxy_Main::p_adapt_enter(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_adapt_enter_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_adapt_enter_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_adapt_enter_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_adapt_enter_void() {
  int epidx = CkRegisterEp("p_adapt_enter()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_enter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_adapt_enter_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_adapt_enter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_adapt_enter_5_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_called();
 */
void CProxy_Main::p_adapt_called(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_adapt_called_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_adapt_called_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_adapt_called_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_adapt_called_void() {
  int epidx = CkRegisterEp("p_adapt_called()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_called_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_adapt_called_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_adapt_called();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_adapt_called_6_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_end();
 */
void CProxy_Main::p_adapt_end(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_adapt_end_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_adapt_end_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_adapt_end_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_adapt_end_void() {
  int epidx = CkRegisterEp("p_adapt_end()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_end_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_adapt_end_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_adapt_end();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_adapt_end_7_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_update();
 */
void CProxy_Main::p_adapt_update(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_adapt_update_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_adapt_update_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_adapt_update_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_adapt_update_void() {
  int epidx = CkRegisterEp("p_adapt_update()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_update_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_adapt_update_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_adapt_update();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_adapt_update_8_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_exit();
 */
void CProxy_Main::p_adapt_exit(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_adapt_exit_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_adapt_exit_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_adapt_exit_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_adapt_exit_void() {
  int epidx = CkRegisterEp("p_adapt_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_adapt_exit_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_adapt_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_adapt_exit_9_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_enter();
 */
void CProxy_Main::p_compute_enter(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_compute_enter_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_compute_enter_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_compute_enter_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_compute_enter_void() {
  int epidx = CkRegisterEp("p_compute_enter()",
      reinterpret_cast<CkCallFnPtr>(_call_p_compute_enter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_compute_enter_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_compute_enter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_compute_enter_10_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_continue();
 */
void CProxy_Main::p_compute_continue(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_compute_continue_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_compute_continue_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_compute_continue_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_compute_continue_void() {
  int epidx = CkRegisterEp("p_compute_continue()",
      reinterpret_cast<CkCallFnPtr>(_call_p_compute_continue_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_compute_continue_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_compute_continue();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_compute_continue_11_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_exit();
 */
void CProxy_Main::p_compute_exit(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_compute_exit_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_compute_exit_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_compute_exit_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_compute_exit_void() {
  int epidx = CkRegisterEp("p_compute_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_compute_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_compute_exit_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_compute_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_compute_exit_12_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_enter();
 */
void CProxy_Main::p_output_enter(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_output_enter_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_output_enter_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_output_enter_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_output_enter_void() {
  int epidx = CkRegisterEp("p_output_enter()",
      reinterpret_cast<CkCallFnPtr>(_call_p_output_enter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_output_enter_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_output_enter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_output_enter_13_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_exit();
 */
void CProxy_Main::p_output_exit(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_output_exit_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_output_exit_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_output_exit_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_output_exit_void() {
  int epidx = CkRegisterEp("p_output_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_output_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_output_exit_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_output_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_output_exit_14_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_enter();
 */
void CProxy_Main::p_stopping_enter(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_stopping_enter_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_stopping_enter_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_stopping_enter_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_stopping_enter_void() {
  int epidx = CkRegisterEp("p_stopping_enter()",
      reinterpret_cast<CkCallFnPtr>(_call_p_stopping_enter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_stopping_enter_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_stopping_enter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_stopping_enter_15_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_balance();
 */
void CProxy_Main::p_stopping_balance(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_stopping_balance_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_stopping_balance_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_stopping_balance_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_stopping_balance_void() {
  int epidx = CkRegisterEp("p_stopping_balance()",
      reinterpret_cast<CkCallFnPtr>(_call_p_stopping_balance_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_stopping_balance_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_stopping_balance();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_stopping_balance_16_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_exit();
 */
void CProxy_Main::p_stopping_exit(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_stopping_exit_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_stopping_exit_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_stopping_exit_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_stopping_exit_void() {
  int epidx = CkRegisterEp("p_stopping_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_stopping_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_stopping_exit_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_stopping_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_stopping_exit_17_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count);
 */
void CProxy_Main::p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_dir, impl_cnt_dir;
  impl_off_dir=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_dir=sizeof(char)*(nd));
  int impl_off_file, impl_cnt_file;
  impl_off_file=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_file=sizeof(char)*(nf));
  int impl_off_line, impl_cnt_line;
  impl_off_line=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_line=sizeof(char)*(nl));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|nd;
    implP|impl_off_dir;
    implP|impl_cnt_dir;
    implP|nf;
    implP|impl_off_file;
    implP|impl_cnt_file;
    implP|nl;
    implP|impl_off_line;
    implP|impl_cnt_line;
    implP|count;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|nd;
    implP|impl_off_dir;
    implP|impl_cnt_dir;
    implP|nf;
    implP|impl_off_file;
    implP|impl_cnt_file;
    implP|nl;
    implP|impl_off_line;
    implP|impl_cnt_line;
    implP|count;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_dir,dir,impl_cnt_dir);
  memcpy(impl_buf+impl_off_file,file,impl_cnt_file);
  memcpy(impl_buf+impl_off_line,line,impl_cnt_line);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_text_file_write_marshall18(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_text_file_write_marshall18(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_text_file_write_marshall18(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_text_file_write_marshall18() {
  int epidx = CkRegisterEp("p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count)",
      reinterpret_cast<CkCallFnPtr>(_call_p_text_file_write_marshall18), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_text_file_write_marshall18);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_text_file_write_marshall18);

  return epidx;
}

void CkIndex_Main::_call_p_text_file_write_marshall18(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> nd;
  implP|nd;
  int impl_off_dir, impl_cnt_dir;
  implP|impl_off_dir;
  implP|impl_cnt_dir;
  PUP::detail::TemporaryObjectHolder<int> nf;
  implP|nf;
  int impl_off_file, impl_cnt_file;
  implP|impl_off_file;
  implP|impl_cnt_file;
  PUP::detail::TemporaryObjectHolder<int> nl;
  implP|nl;
  int impl_off_line, impl_cnt_line;
  implP|impl_off_line;
  implP|impl_cnt_line;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *dir=(char *)(impl_buf+impl_off_dir);
  char *file=(char *)(impl_buf+impl_off_file);
  char *line=(char *)(impl_buf+impl_off_line);
  impl_obj->p_text_file_write(std::move(nd.t), dir, std::move(nf.t), file, std::move(nl.t), line, std::move(count.t));
}
int CkIndex_Main::_callmarshall_p_text_file_write_marshall18(char* impl_buf, void* impl_obj_void) {
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> nd;
  implP|nd;
  int impl_off_dir, impl_cnt_dir;
  implP|impl_off_dir;
  implP|impl_cnt_dir;
  PUP::detail::TemporaryObjectHolder<int> nf;
  implP|nf;
  int impl_off_file, impl_cnt_file;
  implP|impl_off_file;
  implP|impl_cnt_file;
  PUP::detail::TemporaryObjectHolder<int> nl;
  implP|nl;
  int impl_off_line, impl_cnt_line;
  implP|impl_off_line;
  implP|impl_cnt_line;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *dir=(char *)(impl_buf+impl_off_dir);
  char *file=(char *)(impl_buf+impl_off_file);
  char *line=(char *)(impl_buf+impl_off_line);
  impl_obj->p_text_file_write(std::move(nd.t), dir, std::move(nf.t), file, std::move(nl.t), line, std::move(count.t));
  return implP.size();
}
void CkIndex_Main::_marshallmessagepup_p_text_file_write_marshall18(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> nd;
  implP|nd;
  int impl_off_dir, impl_cnt_dir;
  implP|impl_off_dir;
  implP|impl_cnt_dir;
  PUP::detail::TemporaryObjectHolder<int> nf;
  implP|nf;
  int impl_off_file, impl_cnt_file;
  implP|impl_off_file;
  implP|impl_cnt_file;
  PUP::detail::TemporaryObjectHolder<int> nl;
  implP|nl;
  int impl_off_line, impl_cnt_line;
  implP|impl_off_line;
  implP|impl_cnt_line;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *dir=(char *)(impl_buf+impl_off_dir);
  char *file=(char *)(impl_buf+impl_off_file);
  char *line=(char *)(impl_buf+impl_off_line);
  if (implDestP.hasComments()) implDestP.comment("nd");
  implDestP|nd;
  if (implDestP.hasComments()) implDestP.comment("dir");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*dir))<impl_cnt_dir;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|dir[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("nf");
  implDestP|nf;
  if (implDestP.hasComments()) implDestP.comment("file");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*file))<impl_cnt_file;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|file[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("nl");
  implDestP|nl;
  if (implDestP.hasComments()) implDestP.comment("line");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*line))<impl_cnt_line;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|line[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("count");
  implDestP|count;
}
PUPable_def(SINGLE_ARG(Closure_Main::p_text_file_write_18_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_exit();
 */
void CProxy_Main::p_exit(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
    int destPE=CkChareMsgPrep(CkIndex_Main::idx_p_exit_void(), impl_msg, &ckGetChareID());
    if (destPE!=-1) ckDelegatedTo()->ChareSend(ckDelegatedPtr(),CkIndex_Main::idx_p_exit_void(), impl_msg, &ckGetChareID(),destPE);
  } else {
    CkSendMsg(CkIndex_Main::idx_p_exit_void(), impl_msg, &ckGetChareID(),0);
  }
}

// Entry point registration function
int CkIndex_Main::reg_p_exit_void() {
  int epidx = CkRegisterEp("p_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_p_exit_void(void* impl_msg, void* impl_obj_void)
{
  Main* impl_obj = static_cast<Main*>(impl_obj_void);
  impl_obj->p_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Main::p_exit_19_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Main(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_Main::reg_Main_CkMigrateMessage() {
  int epidx = CkRegisterEp("Main(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_Main_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_Main::_call_Main_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  new (impl_obj_void) Main((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Main::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeMainChare);
  CkRegisterBase(__idx, CkIndex_Chare::__idx);
  // REG: Main(CkArgMsg* impl_msg);
  idx_Main_CkArgMsg();
  CkRegisterMainChare(__idx, idx_Main_CkArgMsg());

  // REG: void p_exit(int count_blocks);
  idx_p_exit_marshall2();

  // REG: void p_checkpoint_output(int count, const std::string &dir);
  idx_p_checkpoint_output_marshall3();

  // REG: void p_initial_exit();
  idx_p_initial_exit_void();

  // REG: void p_adapt_enter();
  idx_p_adapt_enter_void();

  // REG: void p_adapt_called();
  idx_p_adapt_called_void();

  // REG: void p_adapt_end();
  idx_p_adapt_end_void();

  // REG: void p_adapt_update();
  idx_p_adapt_update_void();

  // REG: void p_adapt_exit();
  idx_p_adapt_exit_void();

  // REG: void p_compute_enter();
  idx_p_compute_enter_void();

  // REG: void p_compute_continue();
  idx_p_compute_continue_void();

  // REG: void p_compute_exit();
  idx_p_compute_exit_void();

  // REG: void p_output_enter();
  idx_p_output_enter_void();

  // REG: void p_output_exit();
  idx_p_output_exit_void();

  // REG: void p_stopping_enter();
  idx_p_stopping_enter_void();

  // REG: void p_stopping_balance();
  idx_p_stopping_balance_void();

  // REG: void p_stopping_exit();
  idx_p_stopping_exit_void();

  // REG: void p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count);
  idx_p_text_file_write_marshall18();

  // REG: void p_exit();
  idx_p_exit_void();

  // REG: Main(CkMigrateMessage* impl_msg);
  idx_Main_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_Main_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
void _registermain_mesh(void)
{
  static int _done = 0; if(_done) return; _done = 1;
  _registermesh();

  CkRegisterReadonly("proxy_main","CProxy_Main",sizeof(proxy_main),(void *) &proxy_main,__xlater_roPup_proxy_main);

/* REG: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
void p_exit(int count_blocks);
void p_checkpoint_output(int count, const std::string &dir);
void p_initial_exit();
void p_adapt_enter();
void p_adapt_called();
void p_adapt_end();
void p_adapt_update();
void p_adapt_exit();
void p_compute_enter();
void p_compute_continue();
void p_compute_exit();
void p_output_enter();
void p_output_exit();
void p_stopping_enter();
void p_stopping_balance();
void p_stopping_exit();
void p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count);
void p_exit();
Main(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_Main::__register("Main", sizeof(Main));

}
extern "C" void CkRegisterMainModule(void) {
  _registermain_mesh();
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_Main::virtual_pup(PUP::er &p) {
    recursive_pup<Main>(dynamic_cast<Main*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
