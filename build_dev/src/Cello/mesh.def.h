











































































/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_initial_exit_4_closure : public SDAG::Closure {
      

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
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_output_enter_8_closure : public SDAG::Closure {
      

      p_output_enter_8_closure() {
        init();
      }
      p_output_enter_8_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_output_enter_8_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_output_enter_8_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_output_end_10_closure : public SDAG::Closure {
      

      p_output_end_10_closure() {
        init();
      }
      p_output_end_10_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_output_end_10_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_output_end_10_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_output_exit_11_closure : public SDAG::Closure {
      

      p_output_exit_11_closure() {
        init();
      }
      p_output_exit_11_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_output_exit_11_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_output_exit_11_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_output_write_13_closure : public SDAG::Closure {
            int index_output;
            int step;


      p_output_write_13_closure() {
        init();
      }
      p_output_write_13_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return index_output;}
            int & getP1() { return step;}
      void pup(PUP::er& __p) {
        __p | index_output;
        __p | step;
        packClosure(__p);
      }
      virtual ~p_output_write_13_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_output_write_13_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_compute_enter_14_closure : public SDAG::Closure {
      

      p_compute_enter_14_closure() {
        init();
      }
      p_compute_enter_14_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_compute_enter_14_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_compute_enter_14_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_compute_continue_15_closure : public SDAG::Closure {
      

      p_compute_continue_15_closure() {
        init();
      }
      p_compute_continue_15_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_compute_continue_15_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_compute_continue_15_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_compute_exit_17_closure : public SDAG::Closure {
      

      p_compute_exit_17_closure() {
        init();
      }
      p_compute_exit_17_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_compute_exit_17_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_compute_exit_17_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_method_flux_correct_refresh_19_closure : public SDAG::Closure {
      

      p_method_flux_correct_refresh_19_closure() {
        init();
      }
      p_method_flux_correct_refresh_19_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_flux_correct_refresh_19_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_flux_correct_refresh_19_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_method_order_morton_weight_23_closure : public SDAG::Closure {
            int *ic3;
            int weight;
            Index index;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_method_order_morton_weight_23_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_method_order_morton_weight_23_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int *& getP0() { return ic3;}
            int & getP1() { return weight;}
            Index & getP2() { return index;}
      void pup(PUP::er& __p) {
        __p | weight;
        __p | index;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  PUP::detail::TemporaryObjectHolder<int> weight;
  implP|weight;
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
          impl_buf+=CK_ALIGN(implP.size(),16);
          ic3 = (int *)(impl_buf+impl_off_ic3);
        }
      }
      virtual ~p_method_order_morton_weight_23_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_method_order_morton_weight_23_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_method_order_morton_index_24_closure : public SDAG::Closure {
            int index;
            int count;


      p_method_order_morton_index_24_closure() {
        init();
      }
      p_method_order_morton_index_24_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return index;}
            int & getP1() { return count;}
      void pup(PUP::er& __p) {
        __p | index;
        __p | count;
        packClosure(__p);
      }
      virtual ~p_method_order_morton_index_24_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_order_morton_index_24_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_stopping_enter_31_closure : public SDAG::Closure {
      

      p_stopping_enter_31_closure() {
        init();
      }
      p_stopping_enter_31_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_stopping_enter_31_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_stopping_enter_31_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_stopping_load_balance_33_closure : public SDAG::Closure {
      

      p_stopping_load_balance_33_closure() {
        init();
      }
      p_stopping_load_balance_33_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_stopping_load_balance_33_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_stopping_load_balance_33_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_stopping_exit_35_closure : public SDAG::Closure {
      

      p_stopping_exit_35_closure() {
        init();
      }
      p_stopping_exit_35_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_stopping_exit_35_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_stopping_exit_35_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_exit_37_closure : public SDAG::Closure {
      

      p_exit_37_closure() {
        init();
      }
      p_exit_37_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_exit_37_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_exit_37_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_control_sync_count_39_closure : public SDAG::Closure {
            int entry_point;
            int id;
            int count;


      p_control_sync_count_39_closure() {
        init();
      }
      p_control_sync_count_39_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return entry_point;}
            int & getP1() { return id;}
            int & getP2() { return count;}
      void pup(PUP::er& __p) {
        __p | entry_point;
        __p | id;
        __p | count;
        packClosure(__p);
      }
      virtual ~p_control_sync_count_39_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_control_sync_count_39_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_adapt_enter_40_closure : public SDAG::Closure {
      

      p_adapt_enter_40_closure() {
        init();
      }
      p_adapt_enter_40_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_enter_40_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_enter_40_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_adapt_end_42_closure : public SDAG::Closure {
      

      p_adapt_end_42_closure() {
        init();
      }
      p_adapt_end_42_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_end_42_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_end_42_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_adapt_update_43_closure : public SDAG::Closure {
      

      p_adapt_update_43_closure() {
        init();
      }
      p_adapt_update_43_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_update_43_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_update_43_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_adapt_called_45_closure : public SDAG::Closure {
      

      p_adapt_called_45_closure() {
        init();
      }
      p_adapt_called_45_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_called_45_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_called_45_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_adapt_exit_46_closure : public SDAG::Closure {
      

      p_adapt_exit_46_closure() {
        init();
      }
      p_adapt_exit_46_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_exit_46_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_exit_46_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_adapt_delete_47_closure : public SDAG::Closure {
      

      p_adapt_delete_47_closure() {
        init();
      }
      p_adapt_delete_47_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_adapt_delete_47_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_adapt_delete_47_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Block::p_refresh_child_52_closure : public SDAG::Closure {
            int n;
            char *a;
            int *ic3;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_refresh_child_52_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_refresh_child_52_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return n;}
            char *& getP1() { return a;}
            int *& getP2() { return ic3;}
      void pup(PUP::er& __p) {
        __p | n;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_a, impl_cnt_a;
  implP|impl_off_a;
  implP|impl_cnt_a;
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
          impl_buf+=CK_ALIGN(implP.size(),16);
          a = (char *)(impl_buf+impl_off_a);
          ic3 = (int *)(impl_buf+impl_off_ic3);
        }
      }
      virtual ~p_refresh_child_52_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_refresh_child_52_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */











/* DEFS: readonly int MsgCoarsen::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_MsgCoarsen_QColon__QColon_counter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgCoarsen::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgCoarsen::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&MsgCoarsen::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int MsgAdapt::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_MsgAdapt_QColon__QColon_counter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgAdapt::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgAdapt::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&MsgAdapt::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int MsgInitial::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_MsgInitial_QColon__QColon_counter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgInitial::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgInitial::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&MsgInitial::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int MsgOutput::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_MsgOutput_QColon__QColon_counter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgOutput::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgOutput::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&MsgOutput::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int MsgRefine::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_MsgRefine_QColon__QColon_counter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgRefine::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgRefine::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&MsgRefine::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int MsgRefresh::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_MsgRefresh_QColon__QColon_counter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgRefresh::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& MsgRefresh::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&MsgRefresh::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int DataMsg::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_DataMsg_QColon__QColon_counter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& DataMsg::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& DataMsg::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&DataMsg::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int FieldFace::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_FieldFace_QColon__QColon_counter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& FieldFace::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& FieldFace::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&FieldFace::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int ParticleData::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_ParticleData_QColon__QColon_counter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& ParticleData::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& ParticleData::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&ParticleData::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int InitialTrace::id0_[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_InitialTrace_QColon__QColon_id0_(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& InitialTrace::id0_[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      _impl_p|myBuffer;
      if(CkNumNodes() > 1)
        readonlyCreateOnSource(myBuffer);
      PUP::toMem &_impl_p_toMem = *(PUP::toMem *)_impl_pup_er;
      envelope *env = UsrToEnv(_impl_p_toMem.get_orig_pointer());
      CMI_ZC_MSGTYPE(env) = CMK_ZC_BCAST_SEND_MSG;
    }
    if(_impl_p.isUnpacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& InitialTrace::id0_[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
      PUP::fromMem &_impl_p_fromMem = *(PUP::fromMem *)_impl_pup_er;
      char *ptr = _impl_p_fromMem.get_current_pointer();
      PUP::toMem _impl_p_toMem = (PUP::toMem)((void *)ptr);
      envelope *env = UsrToEnv(_impl_p_fromMem.get_orig_pointer());
      CkNcpyBuffer srcBuffer;
      _impl_p|srcBuffer;
      _impl_p_toMem|myBuffer;
      readonlyGet(srcBuffer, myBuffer, (void *)env);
    }
  } else
#endif
  {
  _impl_p(&InitialTrace::id0_[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly double Method::courant_global;
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_Method_QColon__QColon_courant_global(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|Method::courant_global;
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly Config g_config;
 */
extern Config g_config;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_g_config(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|g_config;
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly Parameters g_parameters;
 */
extern Parameters g_parameters;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_g_parameters(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|g_parameters;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(BoundaryPeriodic)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(BoundaryValue)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(ColormapRGB)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Config)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Factory)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Initial)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(InitialTrace)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(InitialValue)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Io)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(IoBlock)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(IoFieldData)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(IoHierarchy)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(IoParticleData)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(IoSimulation)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(ItIndexList)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(ItIndexRange)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MaskExpr)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MaskPng)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MethodCloseFiles)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MethodDebug)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MethodFluxCorrect)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MethodNull)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MethodOrderMorton)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MethodOutput)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MethodRefresh)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(MethodTrace)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(ObjectSphere)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(OutputCheckpoint)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(OutputData)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(OutputImage)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Physics)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Problem)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(ProlongInject)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(ProlongLinear)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Refine)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(RefineDensity)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(RefineMask)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(RefineParticleCount)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(RefineShear)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(RefineSlope)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Refresh)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(RestrictLinear)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(ScheduleInterval)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(ScheduleList)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(SolverNull)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Stopping)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(Units)
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: message FieldMsg{
char a[];
}
;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_FieldMsg::operator new(size_t s){
  return FieldMsg::alloc(__idx, s, 0, 0, GroupDepNum{});
}
void *CMessage_FieldMsg::operator new(size_t s, int* sz){
  return FieldMsg::alloc(__idx, s, sz, 0, GroupDepNum{});
}
void *CMessage_FieldMsg::operator new(size_t s, int* sz,const int pb){
  return FieldMsg::alloc(__idx, s, sz, pb, GroupDepNum{});
}
void *CMessage_FieldMsg::operator new(size_t s, int* sz,const int pb, const GroupDepNum groupDepNum){
  return FieldMsg::alloc(__idx, s, sz, pb, groupDepNum);
}
void *CMessage_FieldMsg::operator new(size_t s, int sz0) {
  int sizes[1];
  sizes[0] = sz0;
  return FieldMsg::alloc(__idx, s, sizes, 0, GroupDepNum{});
}
void *CMessage_FieldMsg::operator new(size_t s, int sz0, const int p) {
  int sizes[1];
  sizes[0] = sz0;
  return FieldMsg::alloc(__idx, s, sizes, p, GroupDepNum{});
}
void *CMessage_FieldMsg::operator new(size_t s, int sz0, const int p, const GroupDepNum groupDepNum) {
  int sizes[1];
  sizes[0] = sz0;
  return FieldMsg::alloc(__idx, s, sizes, p, groupDepNum);
}
void* CMessage_FieldMsg::alloc(int msgnum, size_t sz, int *sizes, int pb, GroupDepNum groupDepNum) {
  CkpvAccess(_offsets)[0] = ALIGN_DEFAULT(sz);
  if(sizes==0)
    CkpvAccess(_offsets)[1] = CkpvAccess(_offsets)[0];
  else
    CkpvAccess(_offsets)[1] = CkpvAccess(_offsets)[0] + ALIGN_DEFAULT(sizeof(char)*sizes[0]);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[1], pb, groupDepNum);
}
CMessage_FieldMsg::CMessage_FieldMsg() {
FieldMsg *newmsg = (FieldMsg *)this;
  newmsg->a = (char *) ((char *)newmsg + CkpvAccess(_offsets)[0]);
}
void CMessage_FieldMsg::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_FieldMsg::pack(FieldMsg *msg) {
  msg->a = (char *) ((char *)msg->a - (char *)msg);
  return (void *) msg;
}
FieldMsg* CMessage_FieldMsg::unpack(void* buf) {
  FieldMsg *msg = (FieldMsg *) buf;
  msg->a = (char *) ((size_t)msg->a + (char *)msg);
  return msg;
}
int CMessage_FieldMsg::__idx=0;
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: message MsgAdapt;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_MsgAdapt::operator new(size_t s){
  return MsgAdapt::alloc(__idx, s, 0, 0, GroupDepNum{});
}
void *CMessage_MsgAdapt::operator new(size_t s, int* sz){
  return MsgAdapt::alloc(__idx, s, sz, 0, GroupDepNum{});
}
void *CMessage_MsgAdapt::operator new(size_t s, int* sz,const int pb){
  return MsgAdapt::alloc(__idx, s, sz, pb, GroupDepNum{});
}
void *CMessage_MsgAdapt::operator new(size_t s, int* sz,const int pb, const GroupDepNum groupDepNum){
  return MsgAdapt::alloc(__idx, s, sz, pb, groupDepNum);
}
void *CMessage_MsgAdapt::operator new(size_t s, const int p) {
  return MsgAdapt::alloc(__idx, s, 0, p, GroupDepNum{});
}
void *CMessage_MsgAdapt::operator new(size_t s, const int p, const GroupDepNum groupDepNum) {
  return MsgAdapt::alloc(__idx, s, 0, p, groupDepNum);
}
void* CMessage_MsgAdapt::alloc(int msgnum, size_t sz, int *sizes, int pb, GroupDepNum groupDepNum) {
  CkpvAccess(_offsets)[0] = ALIGN_DEFAULT(sz);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[0], pb, groupDepNum);
}
CMessage_MsgAdapt::CMessage_MsgAdapt() {
MsgAdapt *newmsg = (MsgAdapt *)this;
}
void CMessage_MsgAdapt::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_MsgAdapt::pack(MsgAdapt *msg) {
  return (void *) msg;
}
MsgAdapt* CMessage_MsgAdapt::unpack(void* buf) {
  MsgAdapt *msg = (MsgAdapt *) buf;
  return msg;
}
int CMessage_MsgAdapt::__idx=0;
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: message MsgCoarsen;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_MsgCoarsen::operator new(size_t s){
  return MsgCoarsen::alloc(__idx, s, 0, 0, GroupDepNum{});
}
void *CMessage_MsgCoarsen::operator new(size_t s, int* sz){
  return MsgCoarsen::alloc(__idx, s, sz, 0, GroupDepNum{});
}
void *CMessage_MsgCoarsen::operator new(size_t s, int* sz,const int pb){
  return MsgCoarsen::alloc(__idx, s, sz, pb, GroupDepNum{});
}
void *CMessage_MsgCoarsen::operator new(size_t s, int* sz,const int pb, const GroupDepNum groupDepNum){
  return MsgCoarsen::alloc(__idx, s, sz, pb, groupDepNum);
}
void *CMessage_MsgCoarsen::operator new(size_t s, const int p) {
  return MsgCoarsen::alloc(__idx, s, 0, p, GroupDepNum{});
}
void *CMessage_MsgCoarsen::operator new(size_t s, const int p, const GroupDepNum groupDepNum) {
  return MsgCoarsen::alloc(__idx, s, 0, p, groupDepNum);
}
void* CMessage_MsgCoarsen::alloc(int msgnum, size_t sz, int *sizes, int pb, GroupDepNum groupDepNum) {
  CkpvAccess(_offsets)[0] = ALIGN_DEFAULT(sz);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[0], pb, groupDepNum);
}
CMessage_MsgCoarsen::CMessage_MsgCoarsen() {
MsgCoarsen *newmsg = (MsgCoarsen *)this;
}
void CMessage_MsgCoarsen::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_MsgCoarsen::pack(MsgCoarsen *msg) {
  return (void *) msg;
}
MsgCoarsen* CMessage_MsgCoarsen::unpack(void* buf) {
  MsgCoarsen *msg = (MsgCoarsen *) buf;
  return msg;
}
int CMessage_MsgCoarsen::__idx=0;
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: message MsgInitial;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_MsgInitial::operator new(size_t s){
  return MsgInitial::alloc(__idx, s, 0, 0, GroupDepNum{});
}
void *CMessage_MsgInitial::operator new(size_t s, int* sz){
  return MsgInitial::alloc(__idx, s, sz, 0, GroupDepNum{});
}
void *CMessage_MsgInitial::operator new(size_t s, int* sz,const int pb){
  return MsgInitial::alloc(__idx, s, sz, pb, GroupDepNum{});
}
void *CMessage_MsgInitial::operator new(size_t s, int* sz,const int pb, const GroupDepNum groupDepNum){
  return MsgInitial::alloc(__idx, s, sz, pb, groupDepNum);
}
void *CMessage_MsgInitial::operator new(size_t s, const int p) {
  return MsgInitial::alloc(__idx, s, 0, p, GroupDepNum{});
}
void *CMessage_MsgInitial::operator new(size_t s, const int p, const GroupDepNum groupDepNum) {
  return MsgInitial::alloc(__idx, s, 0, p, groupDepNum);
}
void* CMessage_MsgInitial::alloc(int msgnum, size_t sz, int *sizes, int pb, GroupDepNum groupDepNum) {
  CkpvAccess(_offsets)[0] = ALIGN_DEFAULT(sz);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[0], pb, groupDepNum);
}
CMessage_MsgInitial::CMessage_MsgInitial() {
MsgInitial *newmsg = (MsgInitial *)this;
}
void CMessage_MsgInitial::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_MsgInitial::pack(MsgInitial *msg) {
  return (void *) msg;
}
MsgInitial* CMessage_MsgInitial::unpack(void* buf) {
  MsgInitial *msg = (MsgInitial *) buf;
  return msg;
}
int CMessage_MsgInitial::__idx=0;
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: message MsgOutput;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_MsgOutput::operator new(size_t s){
  return MsgOutput::alloc(__idx, s, 0, 0, GroupDepNum{});
}
void *CMessage_MsgOutput::operator new(size_t s, int* sz){
  return MsgOutput::alloc(__idx, s, sz, 0, GroupDepNum{});
}
void *CMessage_MsgOutput::operator new(size_t s, int* sz,const int pb){
  return MsgOutput::alloc(__idx, s, sz, pb, GroupDepNum{});
}
void *CMessage_MsgOutput::operator new(size_t s, int* sz,const int pb, const GroupDepNum groupDepNum){
  return MsgOutput::alloc(__idx, s, sz, pb, groupDepNum);
}
void *CMessage_MsgOutput::operator new(size_t s, const int p) {
  return MsgOutput::alloc(__idx, s, 0, p, GroupDepNum{});
}
void *CMessage_MsgOutput::operator new(size_t s, const int p, const GroupDepNum groupDepNum) {
  return MsgOutput::alloc(__idx, s, 0, p, groupDepNum);
}
void* CMessage_MsgOutput::alloc(int msgnum, size_t sz, int *sizes, int pb, GroupDepNum groupDepNum) {
  CkpvAccess(_offsets)[0] = ALIGN_DEFAULT(sz);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[0], pb, groupDepNum);
}
CMessage_MsgOutput::CMessage_MsgOutput() {
MsgOutput *newmsg = (MsgOutput *)this;
}
void CMessage_MsgOutput::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_MsgOutput::pack(MsgOutput *msg) {
  return (void *) msg;
}
MsgOutput* CMessage_MsgOutput::unpack(void* buf) {
  MsgOutput *msg = (MsgOutput *) buf;
  return msg;
}
int CMessage_MsgOutput::__idx=0;
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: message MsgRefine;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_MsgRefine::operator new(size_t s){
  return MsgRefine::alloc(__idx, s, 0, 0, GroupDepNum{});
}
void *CMessage_MsgRefine::operator new(size_t s, int* sz){
  return MsgRefine::alloc(__idx, s, sz, 0, GroupDepNum{});
}
void *CMessage_MsgRefine::operator new(size_t s, int* sz,const int pb){
  return MsgRefine::alloc(__idx, s, sz, pb, GroupDepNum{});
}
void *CMessage_MsgRefine::operator new(size_t s, int* sz,const int pb, const GroupDepNum groupDepNum){
  return MsgRefine::alloc(__idx, s, sz, pb, groupDepNum);
}
void *CMessage_MsgRefine::operator new(size_t s, const int p) {
  return MsgRefine::alloc(__idx, s, 0, p, GroupDepNum{});
}
void *CMessage_MsgRefine::operator new(size_t s, const int p, const GroupDepNum groupDepNum) {
  return MsgRefine::alloc(__idx, s, 0, p, groupDepNum);
}
void* CMessage_MsgRefine::alloc(int msgnum, size_t sz, int *sizes, int pb, GroupDepNum groupDepNum) {
  CkpvAccess(_offsets)[0] = ALIGN_DEFAULT(sz);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[0], pb, groupDepNum);
}
CMessage_MsgRefine::CMessage_MsgRefine() {
MsgRefine *newmsg = (MsgRefine *)this;
}
void CMessage_MsgRefine::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_MsgRefine::pack(MsgRefine *msg) {
  return (void *) msg;
}
MsgRefine* CMessage_MsgRefine::unpack(void* buf) {
  MsgRefine *msg = (MsgRefine *) buf;
  return msg;
}
int CMessage_MsgRefine::__idx=0;
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: message MsgRefresh;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_MsgRefresh::operator new(size_t s){
  return MsgRefresh::alloc(__idx, s, 0, 0, GroupDepNum{});
}
void *CMessage_MsgRefresh::operator new(size_t s, int* sz){
  return MsgRefresh::alloc(__idx, s, sz, 0, GroupDepNum{});
}
void *CMessage_MsgRefresh::operator new(size_t s, int* sz,const int pb){
  return MsgRefresh::alloc(__idx, s, sz, pb, GroupDepNum{});
}
void *CMessage_MsgRefresh::operator new(size_t s, int* sz,const int pb, const GroupDepNum groupDepNum){
  return MsgRefresh::alloc(__idx, s, sz, pb, groupDepNum);
}
void *CMessage_MsgRefresh::operator new(size_t s, const int p) {
  return MsgRefresh::alloc(__idx, s, 0, p, GroupDepNum{});
}
void *CMessage_MsgRefresh::operator new(size_t s, const int p, const GroupDepNum groupDepNum) {
  return MsgRefresh::alloc(__idx, s, 0, p, groupDepNum);
}
void* CMessage_MsgRefresh::alloc(int msgnum, size_t sz, int *sizes, int pb, GroupDepNum groupDepNum) {
  CkpvAccess(_offsets)[0] = ALIGN_DEFAULT(sz);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[0], pb, groupDepNum);
}
CMessage_MsgRefresh::CMessage_MsgRefresh() {
MsgRefresh *newmsg = (MsgRefresh *)this;
}
void CMessage_MsgRefresh::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_MsgRefresh::pack(MsgRefresh *msg) {
  return (void *) msg;
}
MsgRefresh* CMessage_MsgRefresh::unpack(void* buf) {
  MsgRefresh *msg = (MsgRefresh *) buf;
  return msg;
}
int CMessage_MsgRefresh::__idx=0;
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: array Block: ArrayElement{
Block(const process_type &ip_source, const MsgType &msg_type);
void p_set_msg_refine(MsgRefine* impl_msg);
Block();
void p_initial_exit();
void r_end_initialize(CkReductionMsg* impl_msg);
void r_initial_new_next(CkReductionMsg* impl_msg);
void r_initial_new_continue(CkReductionMsg* impl_msg);
void p_output_enter();
void r_output_enter(CkReductionMsg* impl_msg);
void p_output_end();
void p_output_exit();
void r_output_exit(CkReductionMsg* impl_msg);
void p_output_write(int index_output, int step);
void p_compute_enter();
void p_compute_continue();
void r_compute_continue(CkReductionMsg* impl_msg);
void p_compute_exit();
void r_compute_exit(CkReductionMsg* impl_msg);
void p_method_flux_correct_refresh();
void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
void r_method_order_morton_continue(CkReductionMsg* impl_msg);
void r_method_order_morton_complete(CkReductionMsg* impl_msg);
void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
void p_method_order_morton_index(int index, int count);
void p_method_output_next(MsgOutput* impl_msg);
void p_method_output_write(MsgOutput* impl_msg);
void r_method_output_continue(CkReductionMsg* impl_msg);
void r_method_output_done(CkReductionMsg* impl_msg);
void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
void p_stopping_enter();
void r_stopping_enter(CkReductionMsg* impl_msg);
void p_stopping_load_balance();
void r_stopping_load_balance(CkReductionMsg* impl_msg);
void p_stopping_exit();
void r_stopping_exit(CkReductionMsg* impl_msg);
void p_exit();
void r_exit(CkReductionMsg* impl_msg);
void p_control_sync_count(int entry_point, int id, int count);
void p_adapt_enter();
void r_adapt_enter(CkReductionMsg* impl_msg);
void p_adapt_end();
void p_adapt_update();
void r_adapt_next(CkReductionMsg* impl_msg);
void p_adapt_called();
void p_adapt_exit();
void p_adapt_delete();
void p_adapt_recv_level(MsgAdapt* impl_msg);
void p_adapt_recv_child(MsgCoarsen* impl_msg);
void r_restart_enter(CkReductionMsg* impl_msg);
void p_refresh_recv(MsgRefresh* impl_msg);
void p_refresh_child(int n, const char *a, const int *ic3);
Block(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Block::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CProxySection_Block::contribute(CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, userData, fragSize);
}

void CProxySection_Block::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, userData, fragSize);
}

template <typename T>
void CProxySection_Block::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, userData, fragSize);
}

void CProxySection_Block::contribute(CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, cb, userData, fragSize);
}

void CProxySection_Block::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, cb, userData, fragSize);
}

template <typename T>
void CProxySection_Block::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, cb, userData, fragSize);
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Block(const process_type &ip_source, const MsgType &msg_type);
 */
void CProxyElement_Block::insert(const process_type &ip_source, const MsgType &msg_type, int onPE, const CkEntryOptions *impl_e_opts)
{ 
   //Marshall: const process_type &ip_source, const MsgType &msg_type
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<process_type>::type>::type &)ip_source;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<MsgType>::type>::type &)msg_type;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<process_type>::type>::type &)ip_source;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<MsgType>::type>::type &)msg_type;
  }
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_Block::idx_Block_marshall1(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
void CProxyElement_Block::p_set_msg_refine(MsgRefine* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_set_msg_refine_MsgRefine(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Block();
 */
void CProxyElement_Block::insert(int onPE, const CkEntryOptions *impl_e_opts)
{ 
   void *impl_msg = CkAllocSysMsg(impl_e_opts);
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_Block::idx_Block_void(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_initial_exit();
 */
void CProxyElement_Block::p_initial_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_initial_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_end_initialize(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_end_initialize(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_end_initialize_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_initial_new_next(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_initial_new_next(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_initial_new_next_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_initial_new_continue(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_initial_new_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_initial_new_continue_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_enter();
 */
void CProxyElement_Block::p_output_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_output_enter_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_output_enter(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_output_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_output_enter_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_end();
 */
void CProxyElement_Block::p_output_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_output_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_exit();
 */
void CProxyElement_Block::p_output_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_output_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_output_exit(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_output_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_output_exit_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_write(int index_output, int step);
 */
void CProxyElement_Block::p_output_write(int index_output, int step, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int index_output, int step
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index_output;
    implP|step;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index_output;
    implP|step;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_output_write_marshall13(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_enter();
 */
void CProxyElement_Block::p_compute_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_compute_enter_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_continue();
 */
void CProxyElement_Block::p_compute_continue(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_compute_continue_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_compute_continue(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_compute_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_compute_continue_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_exit();
 */
void CProxyElement_Block::p_compute_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_compute_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_compute_exit(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_compute_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_compute_exit_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_flux_correct_refresh();
 */
void CProxyElement_Block::p_method_flux_correct_refresh(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_flux_correct_refresh_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_flux_correct_sum_fields_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_order_morton_continue(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_method_order_morton_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_order_morton_continue_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_order_morton_complete(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_method_order_morton_complete(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_order_morton_complete_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
 */
void CProxyElement_Block::p_method_order_morton_weight(const int *ic3, int weight, const Index &index, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const int *ic3, int weight, const Index &index
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|weight;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|weight;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_order_morton_weight_marshall23(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_order_morton_index(int index, int count);
 */
void CProxyElement_Block::p_method_order_morton_index(int index, int count, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int index, int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index;
    implP|count;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_order_morton_index_marshall24(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_output_next(MsgOutput* impl_msg);
 */
void CProxyElement_Block::p_method_output_next(MsgOutput* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_output_next_MsgOutput(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_output_write(MsgOutput* impl_msg);
 */
void CProxyElement_Block::p_method_output_write(MsgOutput* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_output_write_MsgOutput(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_output_continue(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_method_output_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_output_continue_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_output_done(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_method_output_done(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_output_done_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_method_debug_sum_fields(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_debug_sum_fields_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_stopping_compute_timestep(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_stopping_compute_timestep_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_enter();
 */
void CProxyElement_Block::p_stopping_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_stopping_enter_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_enter(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_stopping_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_stopping_enter_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_load_balance();
 */
void CProxyElement_Block::p_stopping_load_balance(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_stopping_load_balance_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_load_balance(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_stopping_load_balance(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_stopping_load_balance_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_exit();
 */
void CProxyElement_Block::p_stopping_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_stopping_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_exit(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_stopping_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_stopping_exit_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_exit();
 */
void CProxyElement_Block::p_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_exit(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_exit_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_control_sync_count(int entry_point, int id, int count);
 */
void CProxyElement_Block::p_control_sync_count(int entry_point, int id, int count, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int entry_point, int id, int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|entry_point;
    implP|id;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|entry_point;
    implP|id;
    implP|count;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_control_sync_count_marshall39(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_enter();
 */
void CProxyElement_Block::p_adapt_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_enter_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_adapt_enter(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_adapt_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_adapt_enter_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_end();
 */
void CProxyElement_Block::p_adapt_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_update();
 */
void CProxyElement_Block::p_adapt_update(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_update_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_adapt_next(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_adapt_next(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_adapt_next_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_called();
 */
void CProxyElement_Block::p_adapt_called(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_called_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_exit();
 */
void CProxyElement_Block::p_adapt_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_delete();
 */
void CProxyElement_Block::p_adapt_delete(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_delete_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_recv_level(MsgAdapt* impl_msg);
 */
void CProxyElement_Block::p_adapt_recv_level(MsgAdapt* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_recv_level_MsgAdapt(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_recv_child(MsgCoarsen* impl_msg);
 */
void CProxyElement_Block::p_adapt_recv_child(MsgCoarsen* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_recv_child_MsgCoarsen(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_restart_enter(CkReductionMsg* impl_msg);
 */
void CProxyElement_Block::r_restart_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_restart_enter_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_refresh_recv(MsgRefresh* impl_msg);
 */
void CProxyElement_Block::p_refresh_recv(MsgRefresh* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_refresh_recv_MsgRefresh(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_refresh_child(int n, const char *a, const int *ic3);
 */
void CProxyElement_Block::p_refresh_child(int n, const char *a, const int *ic3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const char *a, const int *ic3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_a, impl_cnt_a;
  impl_off_a=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_a=sizeof(char)*(n));
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_a;
    implP|impl_cnt_a;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_a;
    implP|impl_cnt_a;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_a,a,impl_cnt_a);
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_refresh_child_marshall52(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Block(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Block(const process_type &ip_source, const MsgType &msg_type);
 */
CkArrayID CProxy_Block::ckNew(const process_type &ip_source, const MsgType &msg_type, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const process_type &ip_source, const MsgType &msg_type
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<process_type>::type>::type &)ip_source;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<MsgType>::type>::type &)msg_type;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<process_type>::type>::type &)ip_source;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<MsgType>::type>::type &)msg_type;
  }
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_Block::idx_Block_marshall1(), opts);
  return gId;
}
void CProxy_Block::ckNew(const process_type &ip_source, const MsgType &msg_type, const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const process_type &ip_source, const MsgType &msg_type
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<process_type>::type>::type &)ip_source;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<MsgType>::type>::type &)msg_type;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<process_type>::type>::type &)ip_source;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<MsgType>::type>::type &)msg_type;
  }
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_Block::idx_Block_marshall1(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_Block::reg_Block_marshall1() {
  int epidx = CkRegisterEp("Block(const process_type &ip_source, const MsgType &msg_type)",
      reinterpret_cast<CkCallFnPtr>(_call_Block_marshall1), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_Block_marshall1);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_Block_marshall1);

  return epidx;
}

void CkIndex_Block::_call_Block_marshall1(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const process_type &ip_source, const MsgType &msg_type*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<process_type> ip_source;
  implP|ip_source;
  PUP::detail::TemporaryObjectHolder<MsgType> msg_type;
  implP|msg_type;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) Block(std::move(ip_source.t), std::move(msg_type.t));
}
int CkIndex_Block::_callmarshall_Block_marshall1(char* impl_buf, void* impl_obj_void) {
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const process_type &ip_source, const MsgType &msg_type*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<process_type> ip_source;
  implP|ip_source;
  PUP::detail::TemporaryObjectHolder<MsgType> msg_type;
  implP|msg_type;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) Block(std::move(ip_source.t), std::move(msg_type.t));
  return implP.size();
}
void CkIndex_Block::_marshallmessagepup_Block_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const process_type &ip_source, const MsgType &msg_type*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<process_type> ip_source;
  implP|ip_source;
  PUP::detail::TemporaryObjectHolder<MsgType> msg_type;
  implP|msg_type;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("ip_source");
  implDestP|ip_source;
  if (implDestP.hasComments()) implDestP.comment("msg_type");
  implDestP|msg_type;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
void CProxy_Block::p_set_msg_refine(MsgRefine* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_set_msg_refine_MsgRefine(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_set_msg_refine_MsgRefine() {
  int epidx = CkRegisterEp("p_set_msg_refine(MsgRefine* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_set_msg_refine_MsgRefine), CMessage_MsgRefine::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)MsgRefine::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_p_set_msg_refine_MsgRefine(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_set_msg_refine((MsgRefine*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Block();
 */
CkArrayID CProxy_Block::ckNew(const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_Block::idx_Block_void(), opts);
  return gId;
}
void CProxy_Block::ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_Block::idx_Block_void(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_Block::reg_Block_void() {
  int epidx = CkRegisterEp("Block()",
      reinterpret_cast<CkCallFnPtr>(_call_Block_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_Block_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  new (impl_obj_void) Block();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_initial_exit();
 */
void CProxy_Block::p_initial_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_initial_exit_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_initial_exit_void() {
  int epidx = CkRegisterEp("p_initial_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_initial_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_initial_exit_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_initial_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_initial_exit_4_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_end_initialize(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_end_initialize(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_end_initialize_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_end_initialize_CkReductionMsg() {
  int epidx = CkRegisterEp("r_end_initialize(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_end_initialize_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_end_initialize_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_end_initialize((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_initial_new_next(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_initial_new_next(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_initial_new_next_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_initial_new_next_CkReductionMsg() {
  int epidx = CkRegisterEp("r_initial_new_next(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_initial_new_next_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_initial_new_next_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_initial_new_next((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_initial_new_continue(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_initial_new_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_initial_new_continue_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_initial_new_continue_CkReductionMsg() {
  int epidx = CkRegisterEp("r_initial_new_continue(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_initial_new_continue_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_initial_new_continue_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_initial_new_continue((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_enter();
 */
void CProxy_Block::p_output_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_output_enter_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_output_enter_void() {
  int epidx = CkRegisterEp("p_output_enter()",
      reinterpret_cast<CkCallFnPtr>(_call_p_output_enter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_output_enter_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_output_enter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_output_enter_8_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_output_enter(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_output_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_output_enter_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_output_enter_CkReductionMsg() {
  int epidx = CkRegisterEp("r_output_enter(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_output_enter_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_output_enter_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_output_enter((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_end();
 */
void CProxy_Block::p_output_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_output_end_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_output_end_void() {
  int epidx = CkRegisterEp("p_output_end()",
      reinterpret_cast<CkCallFnPtr>(_call_p_output_end_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_output_end_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_output_end();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_output_end_10_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_exit();
 */
void CProxy_Block::p_output_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_output_exit_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_output_exit_void() {
  int epidx = CkRegisterEp("p_output_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_output_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_output_exit_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_output_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_output_exit_11_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_output_exit(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_output_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_output_exit_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_output_exit_CkReductionMsg() {
  int epidx = CkRegisterEp("r_output_exit(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_output_exit_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_output_exit_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_output_exit((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_write(int index_output, int step);
 */
void CProxy_Block::p_output_write(int index_output, int step, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int index_output, int step
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index_output;
    implP|step;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index_output;
    implP|step;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_output_write_marshall13(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_output_write_marshall13() {
  int epidx = CkRegisterEp("p_output_write(int index_output, int step)",
      reinterpret_cast<CkCallFnPtr>(_call_p_output_write_marshall13), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_output_write_marshall13);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_output_write_marshall13);

  return epidx;
}

void CkIndex_Block::_call_p_output_write_marshall13(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int index_output, int step*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> index_output;
  implP|index_output;
  PUP::detail::TemporaryObjectHolder<int> step;
  implP|step;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_output_write(std::move(index_output.t), std::move(step.t));
}
int CkIndex_Block::_callmarshall_p_output_write_marshall13(char* impl_buf, void* impl_obj_void) {
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int index_output, int step*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> index_output;
  implP|index_output;
  PUP::detail::TemporaryObjectHolder<int> step;
  implP|step;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_output_write(std::move(index_output.t), std::move(step.t));
  return implP.size();
}
void CkIndex_Block::_marshallmessagepup_p_output_write_marshall13(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int index_output, int step*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> index_output;
  implP|index_output;
  PUP::detail::TemporaryObjectHolder<int> step;
  implP|step;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("index_output");
  implDestP|index_output;
  if (implDestP.hasComments()) implDestP.comment("step");
  implDestP|step;
}
PUPable_def(SINGLE_ARG(Closure_Block::p_output_write_13_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_enter();
 */
void CProxy_Block::p_compute_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_compute_enter_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_compute_enter_void() {
  int epidx = CkRegisterEp("p_compute_enter()",
      reinterpret_cast<CkCallFnPtr>(_call_p_compute_enter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_compute_enter_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_compute_enter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_compute_enter_14_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_continue();
 */
void CProxy_Block::p_compute_continue(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_compute_continue_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_compute_continue_void() {
  int epidx = CkRegisterEp("p_compute_continue()",
      reinterpret_cast<CkCallFnPtr>(_call_p_compute_continue_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_compute_continue_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_compute_continue();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_compute_continue_15_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_compute_continue(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_compute_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_compute_continue_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_compute_continue_CkReductionMsg() {
  int epidx = CkRegisterEp("r_compute_continue(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_compute_continue_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_compute_continue_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_compute_continue((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_exit();
 */
void CProxy_Block::p_compute_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_compute_exit_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_compute_exit_void() {
  int epidx = CkRegisterEp("p_compute_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_compute_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_compute_exit_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_compute_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_compute_exit_17_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_compute_exit(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_compute_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_compute_exit_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_compute_exit_CkReductionMsg() {
  int epidx = CkRegisterEp("r_compute_exit(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_compute_exit_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_compute_exit_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_compute_exit((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_flux_correct_refresh();
 */
void CProxy_Block::p_method_flux_correct_refresh(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_method_flux_correct_refresh_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_method_flux_correct_refresh_void() {
  int epidx = CkRegisterEp("p_method_flux_correct_refresh()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_flux_correct_refresh_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_method_flux_correct_refresh_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_method_flux_correct_refresh();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_method_flux_correct_refresh_19_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_method_flux_correct_sum_fields_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_method_flux_correct_sum_fields_CkReductionMsg() {
  int epidx = CkRegisterEp("r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_method_flux_correct_sum_fields_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_method_flux_correct_sum_fields_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_method_flux_correct_sum_fields((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_order_morton_continue(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_method_order_morton_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_method_order_morton_continue_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_method_order_morton_continue_CkReductionMsg() {
  int epidx = CkRegisterEp("r_method_order_morton_continue(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_method_order_morton_continue_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_method_order_morton_continue_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_method_order_morton_continue((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_order_morton_complete(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_method_order_morton_complete(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_method_order_morton_complete_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_method_order_morton_complete_CkReductionMsg() {
  int epidx = CkRegisterEp("r_method_order_morton_complete(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_method_order_morton_complete_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_method_order_morton_complete_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_method_order_morton_complete((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
 */
void CProxy_Block::p_method_order_morton_weight(const int *ic3, int weight, const Index &index, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const int *ic3, int weight, const Index &index
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|weight;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|weight;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_method_order_morton_weight_marshall23(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_method_order_morton_weight_marshall23() {
  int epidx = CkRegisterEp("p_method_order_morton_weight(const int *ic3, int weight, const Index &index)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_order_morton_weight_marshall23), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_method_order_morton_weight_marshall23);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_method_order_morton_weight_marshall23);

  return epidx;
}

void CkIndex_Block::_call_p_method_order_morton_weight_marshall23(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const int *ic3, int weight, const Index &index*/
  PUP::fromMem implP(impl_buf);
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  PUP::detail::TemporaryObjectHolder<int> weight;
  implP|weight;
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  impl_obj->p_method_order_morton_weight(ic3, std::move(weight.t), std::move(index.t));
}
int CkIndex_Block::_callmarshall_p_method_order_morton_weight_marshall23(char* impl_buf, void* impl_obj_void) {
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const int *ic3, int weight, const Index &index*/
  PUP::fromMem implP(impl_buf);
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  PUP::detail::TemporaryObjectHolder<int> weight;
  implP|weight;
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  impl_obj->p_method_order_morton_weight(ic3, std::move(weight.t), std::move(index.t));
  return implP.size();
}
void CkIndex_Block::_marshallmessagepup_p_method_order_morton_weight_marshall23(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const int *ic3, int weight, const Index &index*/
  PUP::fromMem implP(impl_buf);
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  PUP::detail::TemporaryObjectHolder<int> weight;
  implP|weight;
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  if (implDestP.hasComments()) implDestP.comment("ic3");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*ic3))<impl_cnt_ic3;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|ic3[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("weight");
  implDestP|weight;
  if (implDestP.hasComments()) implDestP.comment("index");
  implDestP|index;
}
PUPable_def(SINGLE_ARG(Closure_Block::p_method_order_morton_weight_23_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_order_morton_index(int index, int count);
 */
void CProxy_Block::p_method_order_morton_index(int index, int count, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int index, int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index;
    implP|count;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_method_order_morton_index_marshall24(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_method_order_morton_index_marshall24() {
  int epidx = CkRegisterEp("p_method_order_morton_index(int index, int count)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_order_morton_index_marshall24), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_method_order_morton_index_marshall24);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_method_order_morton_index_marshall24);

  return epidx;
}

void CkIndex_Block::_call_p_method_order_morton_index_marshall24(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int index, int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> index;
  implP|index;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_method_order_morton_index(std::move(index.t), std::move(count.t));
}
int CkIndex_Block::_callmarshall_p_method_order_morton_index_marshall24(char* impl_buf, void* impl_obj_void) {
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int index, int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> index;
  implP|index;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_method_order_morton_index(std::move(index.t), std::move(count.t));
  return implP.size();
}
void CkIndex_Block::_marshallmessagepup_p_method_order_morton_index_marshall24(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int index, int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> index;
  implP|index;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("index");
  implDestP|index;
  if (implDestP.hasComments()) implDestP.comment("count");
  implDestP|count;
}
PUPable_def(SINGLE_ARG(Closure_Block::p_method_order_morton_index_24_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_output_next(MsgOutput* impl_msg);
 */
void CProxy_Block::p_method_output_next(MsgOutput* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_method_output_next_MsgOutput(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_method_output_next_MsgOutput() {
  int epidx = CkRegisterEp("p_method_output_next(MsgOutput* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_output_next_MsgOutput), CMessage_MsgOutput::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)MsgOutput::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_p_method_output_next_MsgOutput(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_method_output_next((MsgOutput*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_output_write(MsgOutput* impl_msg);
 */
void CProxy_Block::p_method_output_write(MsgOutput* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_method_output_write_MsgOutput(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_method_output_write_MsgOutput() {
  int epidx = CkRegisterEp("p_method_output_write(MsgOutput* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_output_write_MsgOutput), CMessage_MsgOutput::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)MsgOutput::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_p_method_output_write_MsgOutput(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_method_output_write((MsgOutput*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_output_continue(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_method_output_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_method_output_continue_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_method_output_continue_CkReductionMsg() {
  int epidx = CkRegisterEp("r_method_output_continue(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_method_output_continue_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_method_output_continue_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_method_output_continue((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_output_done(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_method_output_done(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_method_output_done_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_method_output_done_CkReductionMsg() {
  int epidx = CkRegisterEp("r_method_output_done(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_method_output_done_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_method_output_done_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_method_output_done((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_method_debug_sum_fields(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_method_debug_sum_fields_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_method_debug_sum_fields_CkReductionMsg() {
  int epidx = CkRegisterEp("r_method_debug_sum_fields(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_method_debug_sum_fields_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_method_debug_sum_fields_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_method_debug_sum_fields((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_stopping_compute_timestep(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_stopping_compute_timestep_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_stopping_compute_timestep_CkReductionMsg() {
  int epidx = CkRegisterEp("r_stopping_compute_timestep(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_stopping_compute_timestep_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_stopping_compute_timestep_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_stopping_compute_timestep((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_enter();
 */
void CProxy_Block::p_stopping_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_stopping_enter_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_stopping_enter_void() {
  int epidx = CkRegisterEp("p_stopping_enter()",
      reinterpret_cast<CkCallFnPtr>(_call_p_stopping_enter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_stopping_enter_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_stopping_enter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_stopping_enter_31_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_enter(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_stopping_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_stopping_enter_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_stopping_enter_CkReductionMsg() {
  int epidx = CkRegisterEp("r_stopping_enter(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_stopping_enter_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_stopping_enter_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_stopping_enter((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_load_balance();
 */
void CProxy_Block::p_stopping_load_balance(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_stopping_load_balance_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_stopping_load_balance_void() {
  int epidx = CkRegisterEp("p_stopping_load_balance()",
      reinterpret_cast<CkCallFnPtr>(_call_p_stopping_load_balance_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_stopping_load_balance_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_stopping_load_balance();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_stopping_load_balance_33_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_load_balance(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_stopping_load_balance(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_stopping_load_balance_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_stopping_load_balance_CkReductionMsg() {
  int epidx = CkRegisterEp("r_stopping_load_balance(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_stopping_load_balance_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_stopping_load_balance_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_stopping_load_balance((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_exit();
 */
void CProxy_Block::p_stopping_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_stopping_exit_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_stopping_exit_void() {
  int epidx = CkRegisterEp("p_stopping_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_stopping_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_stopping_exit_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_stopping_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_stopping_exit_35_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_exit(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_stopping_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_stopping_exit_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_stopping_exit_CkReductionMsg() {
  int epidx = CkRegisterEp("r_stopping_exit(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_stopping_exit_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_stopping_exit_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_stopping_exit((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_exit();
 */
void CProxy_Block::p_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_exit_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_exit_void() {
  int epidx = CkRegisterEp("p_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_exit_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_exit_37_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_exit(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_exit_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_exit_CkReductionMsg() {
  int epidx = CkRegisterEp("r_exit(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_exit_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_exit_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_exit((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_control_sync_count(int entry_point, int id, int count);
 */
void CProxy_Block::p_control_sync_count(int entry_point, int id, int count, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int entry_point, int id, int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|entry_point;
    implP|id;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|entry_point;
    implP|id;
    implP|count;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_control_sync_count_marshall39(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_control_sync_count_marshall39() {
  int epidx = CkRegisterEp("p_control_sync_count(int entry_point, int id, int count)",
      reinterpret_cast<CkCallFnPtr>(_call_p_control_sync_count_marshall39), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_control_sync_count_marshall39);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_control_sync_count_marshall39);

  return epidx;
}

void CkIndex_Block::_call_p_control_sync_count_marshall39(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int entry_point, int id, int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> entry_point;
  implP|entry_point;
  PUP::detail::TemporaryObjectHolder<int> id;
  implP|id;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_control_sync_count(std::move(entry_point.t), std::move(id.t), std::move(count.t));
}
int CkIndex_Block::_callmarshall_p_control_sync_count_marshall39(char* impl_buf, void* impl_obj_void) {
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int entry_point, int id, int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> entry_point;
  implP|entry_point;
  PUP::detail::TemporaryObjectHolder<int> id;
  implP|id;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_control_sync_count(std::move(entry_point.t), std::move(id.t), std::move(count.t));
  return implP.size();
}
void CkIndex_Block::_marshallmessagepup_p_control_sync_count_marshall39(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int entry_point, int id, int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> entry_point;
  implP|entry_point;
  PUP::detail::TemporaryObjectHolder<int> id;
  implP|id;
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("entry_point");
  implDestP|entry_point;
  if (implDestP.hasComments()) implDestP.comment("id");
  implDestP|id;
  if (implDestP.hasComments()) implDestP.comment("count");
  implDestP|count;
}
PUPable_def(SINGLE_ARG(Closure_Block::p_control_sync_count_39_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_enter();
 */
void CProxy_Block::p_adapt_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_adapt_enter_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_adapt_enter_void() {
  int epidx = CkRegisterEp("p_adapt_enter()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_enter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_adapt_enter_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_adapt_enter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_adapt_enter_40_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_adapt_enter(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_adapt_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_adapt_enter_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_adapt_enter_CkReductionMsg() {
  int epidx = CkRegisterEp("r_adapt_enter(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_adapt_enter_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_adapt_enter_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_adapt_enter((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_end();
 */
void CProxy_Block::p_adapt_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_adapt_end_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_adapt_end_void() {
  int epidx = CkRegisterEp("p_adapt_end()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_end_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_adapt_end_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_adapt_end();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_adapt_end_42_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_update();
 */
void CProxy_Block::p_adapt_update(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_adapt_update_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_adapt_update_void() {
  int epidx = CkRegisterEp("p_adapt_update()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_update_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_adapt_update_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_adapt_update();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_adapt_update_43_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_adapt_next(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_adapt_next(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_adapt_next_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_adapt_next_CkReductionMsg() {
  int epidx = CkRegisterEp("r_adapt_next(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_adapt_next_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_adapt_next_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_adapt_next((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_called();
 */
void CProxy_Block::p_adapt_called(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_adapt_called_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_adapt_called_void() {
  int epidx = CkRegisterEp("p_adapt_called()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_called_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_adapt_called_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_adapt_called();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_adapt_called_45_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_exit();
 */
void CProxy_Block::p_adapt_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_adapt_exit_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_adapt_exit_void() {
  int epidx = CkRegisterEp("p_adapt_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_adapt_exit_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_adapt_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_adapt_exit_46_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_delete();
 */
void CProxy_Block::p_adapt_delete(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_adapt_delete_void(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_adapt_delete_void() {
  int epidx = CkRegisterEp("p_adapt_delete()",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_delete_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_p_adapt_delete_void(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_adapt_delete();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_adapt_delete_47_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_recv_level(MsgAdapt* impl_msg);
 */
void CProxy_Block::p_adapt_recv_level(MsgAdapt* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_adapt_recv_level_MsgAdapt(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_adapt_recv_level_MsgAdapt() {
  int epidx = CkRegisterEp("p_adapt_recv_level(MsgAdapt* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_recv_level_MsgAdapt), CMessage_MsgAdapt::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)MsgAdapt::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_p_adapt_recv_level_MsgAdapt(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_adapt_recv_level((MsgAdapt*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_recv_child(MsgCoarsen* impl_msg);
 */
void CProxy_Block::p_adapt_recv_child(MsgCoarsen* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_adapt_recv_child_MsgCoarsen(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_adapt_recv_child_MsgCoarsen() {
  int epidx = CkRegisterEp("p_adapt_recv_child(MsgCoarsen* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_adapt_recv_child_MsgCoarsen), CMessage_MsgCoarsen::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)MsgCoarsen::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_p_adapt_recv_child_MsgCoarsen(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_adapt_recv_child((MsgCoarsen*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_restart_enter(CkReductionMsg* impl_msg);
 */
void CProxy_Block::r_restart_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_r_restart_enter_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_Block::reg_r_restart_enter_CkReductionMsg() {
  int epidx = CkRegisterEp("r_restart_enter(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_restart_enter_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_r_restart_enter_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->r_restart_enter((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_refresh_recv(MsgRefresh* impl_msg);
 */
void CProxy_Block::p_refresh_recv(MsgRefresh* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_refresh_recv_MsgRefresh(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_refresh_recv_MsgRefresh() {
  int epidx = CkRegisterEp("p_refresh_recv(MsgRefresh* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_refresh_recv_MsgRefresh), CMessage_MsgRefresh::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)MsgRefresh::ckDebugPup);
  return epidx;
}

void CkIndex_Block::_call_p_refresh_recv_MsgRefresh(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  impl_obj->p_refresh_recv((MsgRefresh*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_refresh_child(int n, const char *a, const int *ic3);
 */
void CProxy_Block::p_refresh_child(int n, const char *a, const int *ic3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const char *a, const int *ic3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_a, impl_cnt_a;
  impl_off_a=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_a=sizeof(char)*(n));
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_a;
    implP|impl_cnt_a;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_a;
    implP|impl_cnt_a;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_a,a,impl_cnt_a);
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_Block::idx_p_refresh_child_marshall52(),0);
}

// Entry point registration function
int CkIndex_Block::reg_p_refresh_child_marshall52() {
  int epidx = CkRegisterEp("p_refresh_child(int n, const char *a, const int *ic3)",
      reinterpret_cast<CkCallFnPtr>(_call_p_refresh_child_marshall52), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_refresh_child_marshall52);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_refresh_child_marshall52);

  return epidx;
}

void CkIndex_Block::_call_p_refresh_child_marshall52(void* impl_msg, void* impl_obj_void)
{
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const char *a, const int *ic3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_a, impl_cnt_a;
  implP|impl_off_a;
  implP|impl_cnt_a;
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *a=(char *)(impl_buf+impl_off_a);
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  impl_obj->p_refresh_child(std::move(n.t), a, ic3);
}
int CkIndex_Block::_callmarshall_p_refresh_child_marshall52(char* impl_buf, void* impl_obj_void) {
  Block* impl_obj = static_cast<Block*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int n, const char *a, const int *ic3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_a, impl_cnt_a;
  implP|impl_off_a;
  implP|impl_cnt_a;
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *a=(char *)(impl_buf+impl_off_a);
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  impl_obj->p_refresh_child(std::move(n.t), a, ic3);
  return implP.size();
}
void CkIndex_Block::_marshallmessagepup_p_refresh_child_marshall52(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const char *a, const int *ic3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_a, impl_cnt_a;
  implP|impl_off_a;
  implP|impl_cnt_a;
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *a=(char *)(impl_buf+impl_off_a);
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  if (implDestP.hasComments()) implDestP.comment("n");
  implDestP|n;
  if (implDestP.hasComments()) implDestP.comment("a");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*a))<impl_cnt_a;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|a[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("ic3");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*ic3))<impl_cnt_ic3;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|ic3[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Block::p_refresh_child_52_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Block(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_Block::reg_Block_CkMigrateMessage() {
  int epidx = CkRegisterEp("Block(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_Block_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_Block::_call_Block_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  call_migration_constructor<Block> c = impl_obj_void;
  c((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Block(const process_type &ip_source, const MsgType &msg_type);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
void CProxySection_Block::p_set_msg_refine(MsgRefine* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_set_msg_refine_MsgRefine(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Block();
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_initial_exit();
 */
void CProxySection_Block::p_initial_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_initial_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_end_initialize(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_end_initialize(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_end_initialize_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_initial_new_next(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_initial_new_next(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_initial_new_next_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_initial_new_continue(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_initial_new_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_initial_new_continue_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_enter();
 */
void CProxySection_Block::p_output_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_output_enter_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_output_enter(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_output_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_output_enter_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_end();
 */
void CProxySection_Block::p_output_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_output_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_exit();
 */
void CProxySection_Block::p_output_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_output_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_output_exit(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_output_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_output_exit_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_write(int index_output, int step);
 */
void CProxySection_Block::p_output_write(int index_output, int step, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int index_output, int step
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index_output;
    implP|step;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index_output;
    implP|step;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_output_write_marshall13(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_enter();
 */
void CProxySection_Block::p_compute_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_compute_enter_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_continue();
 */
void CProxySection_Block::p_compute_continue(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_compute_continue_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_compute_continue(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_compute_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_compute_continue_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_compute_exit();
 */
void CProxySection_Block::p_compute_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_compute_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_compute_exit(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_compute_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_compute_exit_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_flux_correct_refresh();
 */
void CProxySection_Block::p_method_flux_correct_refresh(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_flux_correct_refresh_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_flux_correct_sum_fields_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_order_morton_continue(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_method_order_morton_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_order_morton_continue_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_order_morton_complete(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_method_order_morton_complete(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_order_morton_complete_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
 */
void CProxySection_Block::p_method_order_morton_weight(const int *ic3, int weight, const Index &index, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const int *ic3, int weight, const Index &index
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|weight;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|weight;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_order_morton_weight_marshall23(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_order_morton_index(int index, int count);
 */
void CProxySection_Block::p_method_order_morton_index(int index, int count, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int index, int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index;
    implP|count;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_order_morton_index_marshall24(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_output_next(MsgOutput* impl_msg);
 */
void CProxySection_Block::p_method_output_next(MsgOutput* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_output_next_MsgOutput(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_output_write(MsgOutput* impl_msg);
 */
void CProxySection_Block::p_method_output_write(MsgOutput* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_method_output_write_MsgOutput(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_output_continue(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_method_output_continue(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_output_continue_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_output_done(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_method_output_done(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_output_done_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_method_debug_sum_fields(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_method_debug_sum_fields_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_stopping_compute_timestep(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_stopping_compute_timestep_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_enter();
 */
void CProxySection_Block::p_stopping_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_stopping_enter_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_enter(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_stopping_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_stopping_enter_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_load_balance();
 */
void CProxySection_Block::p_stopping_load_balance(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_stopping_load_balance_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_load_balance(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_stopping_load_balance(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_stopping_load_balance_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_stopping_exit();
 */
void CProxySection_Block::p_stopping_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_stopping_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_stopping_exit(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_stopping_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_stopping_exit_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_exit();
 */
void CProxySection_Block::p_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_exit(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_exit(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_exit_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_control_sync_count(int entry_point, int id, int count);
 */
void CProxySection_Block::p_control_sync_count(int entry_point, int id, int count, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int entry_point, int id, int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|entry_point;
    implP|id;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|entry_point;
    implP|id;
    implP|count;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_control_sync_count_marshall39(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_enter();
 */
void CProxySection_Block::p_adapt_enter(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_enter_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_adapt_enter(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_adapt_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_adapt_enter_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_end();
 */
void CProxySection_Block::p_adapt_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_update();
 */
void CProxySection_Block::p_adapt_update(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_update_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_adapt_next(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_adapt_next(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_adapt_next_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_called();
 */
void CProxySection_Block::p_adapt_called(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_called_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_exit();
 */
void CProxySection_Block::p_adapt_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_delete();
 */
void CProxySection_Block::p_adapt_delete(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_delete_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_recv_level(MsgAdapt* impl_msg);
 */
void CProxySection_Block::p_adapt_recv_level(MsgAdapt* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_recv_level_MsgAdapt(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_adapt_recv_child(MsgCoarsen* impl_msg);
 */
void CProxySection_Block::p_adapt_recv_child(MsgCoarsen* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_adapt_recv_child_MsgCoarsen(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_restart_enter(CkReductionMsg* impl_msg);
 */
void CProxySection_Block::r_restart_enter(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_r_restart_enter_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_refresh_recv(MsgRefresh* impl_msg);
 */
void CProxySection_Block::p_refresh_recv(MsgRefresh* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_refresh_recv_MsgRefresh(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_refresh_child(int n, const char *a, const int *ic3);
 */
void CProxySection_Block::p_refresh_child(int n, const char *a, const int *ic3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const char *a, const int *ic3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_a, impl_cnt_a;
  impl_off_a=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_a=sizeof(char)*(n));
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_a;
    implP|impl_cnt_a;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_a;
    implP|impl_cnt_a;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_a,a,impl_cnt_a);
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_Block::idx_p_refresh_child_marshall52(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Block(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Block::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterArrayDimensions(__idx, -1);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
  // REG: Block(const process_type &ip_source, const MsgType &msg_type);
  idx_Block_marshall1();

  // REG: void p_set_msg_refine(MsgRefine* impl_msg);
  idx_p_set_msg_refine_MsgRefine();

  // REG: Block();
  idx_Block_void();
  CkRegisterDefaultCtor(__idx, idx_Block_void());

  // REG: void p_initial_exit();
  idx_p_initial_exit_void();

  // REG: void r_end_initialize(CkReductionMsg* impl_msg);
  idx_r_end_initialize_CkReductionMsg();

  // REG: void r_initial_new_next(CkReductionMsg* impl_msg);
  idx_r_initial_new_next_CkReductionMsg();

  // REG: void r_initial_new_continue(CkReductionMsg* impl_msg);
  idx_r_initial_new_continue_CkReductionMsg();

  // REG: void p_output_enter();
  idx_p_output_enter_void();

  // REG: void r_output_enter(CkReductionMsg* impl_msg);
  idx_r_output_enter_CkReductionMsg();

  // REG: void p_output_end();
  idx_p_output_end_void();

  // REG: void p_output_exit();
  idx_p_output_exit_void();

  // REG: void r_output_exit(CkReductionMsg* impl_msg);
  idx_r_output_exit_CkReductionMsg();

  // REG: void p_output_write(int index_output, int step);
  idx_p_output_write_marshall13();

  // REG: void p_compute_enter();
  idx_p_compute_enter_void();

  // REG: void p_compute_continue();
  idx_p_compute_continue_void();

  // REG: void r_compute_continue(CkReductionMsg* impl_msg);
  idx_r_compute_continue_CkReductionMsg();

  // REG: void p_compute_exit();
  idx_p_compute_exit_void();

  // REG: void r_compute_exit(CkReductionMsg* impl_msg);
  idx_r_compute_exit_CkReductionMsg();

  // REG: void p_method_flux_correct_refresh();
  idx_p_method_flux_correct_refresh_void();

  // REG: void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
  idx_r_method_flux_correct_sum_fields_CkReductionMsg();

  // REG: void r_method_order_morton_continue(CkReductionMsg* impl_msg);
  idx_r_method_order_morton_continue_CkReductionMsg();

  // REG: void r_method_order_morton_complete(CkReductionMsg* impl_msg);
  idx_r_method_order_morton_complete_CkReductionMsg();

  // REG: void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
  idx_p_method_order_morton_weight_marshall23();

  // REG: void p_method_order_morton_index(int index, int count);
  idx_p_method_order_morton_index_marshall24();

  // REG: void p_method_output_next(MsgOutput* impl_msg);
  idx_p_method_output_next_MsgOutput();

  // REG: void p_method_output_write(MsgOutput* impl_msg);
  idx_p_method_output_write_MsgOutput();

  // REG: void r_method_output_continue(CkReductionMsg* impl_msg);
  idx_r_method_output_continue_CkReductionMsg();

  // REG: void r_method_output_done(CkReductionMsg* impl_msg);
  idx_r_method_output_done_CkReductionMsg();

  // REG: void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
  idx_r_method_debug_sum_fields_CkReductionMsg();

  // REG: void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
  idx_r_stopping_compute_timestep_CkReductionMsg();

  // REG: void p_stopping_enter();
  idx_p_stopping_enter_void();

  // REG: void r_stopping_enter(CkReductionMsg* impl_msg);
  idx_r_stopping_enter_CkReductionMsg();

  // REG: void p_stopping_load_balance();
  idx_p_stopping_load_balance_void();

  // REG: void r_stopping_load_balance(CkReductionMsg* impl_msg);
  idx_r_stopping_load_balance_CkReductionMsg();

  // REG: void p_stopping_exit();
  idx_p_stopping_exit_void();

  // REG: void r_stopping_exit(CkReductionMsg* impl_msg);
  idx_r_stopping_exit_CkReductionMsg();

  // REG: void p_exit();
  idx_p_exit_void();

  // REG: void r_exit(CkReductionMsg* impl_msg);
  idx_r_exit_CkReductionMsg();

  // REG: void p_control_sync_count(int entry_point, int id, int count);
  idx_p_control_sync_count_marshall39();

  // REG: void p_adapt_enter();
  idx_p_adapt_enter_void();

  // REG: void r_adapt_enter(CkReductionMsg* impl_msg);
  idx_r_adapt_enter_CkReductionMsg();

  // REG: void p_adapt_end();
  idx_p_adapt_end_void();

  // REG: void p_adapt_update();
  idx_p_adapt_update_void();

  // REG: void r_adapt_next(CkReductionMsg* impl_msg);
  idx_r_adapt_next_CkReductionMsg();

  // REG: void p_adapt_called();
  idx_p_adapt_called_void();

  // REG: void p_adapt_exit();
  idx_p_adapt_exit_void();

  // REG: void p_adapt_delete();
  idx_p_adapt_delete_void();

  // REG: void p_adapt_recv_level(MsgAdapt* impl_msg);
  idx_p_adapt_recv_level_MsgAdapt();

  // REG: void p_adapt_recv_child(MsgCoarsen* impl_msg);
  idx_p_adapt_recv_child_MsgCoarsen();

  // REG: void r_restart_enter(CkReductionMsg* impl_msg);
  idx_r_restart_enter_CkReductionMsg();

  // REG: void p_refresh_recv(MsgRefresh* impl_msg);
  idx_p_refresh_recv_MsgRefresh();

  // REG: void p_refresh_child(int n, const char *a, const int *ic3);
  idx_p_refresh_child_marshall52();

  // REG: Block(CkMigrateMessage* impl_msg);
  idx_Block_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_Block_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: array IoReader: ArrayElement{
IoReader();
IoReader(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_IoReader::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CProxySection_IoReader::contribute(CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, userData, fragSize);
}

void CProxySection_IoReader::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, userData, fragSize);
}

template <typename T>
void CProxySection_IoReader::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, userData, fragSize);
}

void CProxySection_IoReader::contribute(CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, cb, userData, fragSize);
}

void CProxySection_IoReader::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, cb, userData, fragSize);
}

template <typename T>
void CProxySection_IoReader::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, cb, userData, fragSize);
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoReader();
 */
void CProxyElement_IoReader::insert(int onPE, const CkEntryOptions *impl_e_opts)
{ 
   void *impl_msg = CkAllocSysMsg(impl_e_opts);
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_IoReader::idx_IoReader_void(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoReader(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoReader();
 */
CkArrayID CProxy_IoReader::ckNew(const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoReader::idx_IoReader_void(), opts);
  return gId;
}
void CProxy_IoReader::ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoReader::idx_IoReader_void(), _ck_array_creation_cb, opts, impl_msg);
}
CkArrayID CProxy_IoReader::ckNew(const int s1, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoReader::idx_IoReader_void(), opts);
  return gId;
}
void CProxy_IoReader::ckNew(const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoReader::idx_IoReader_void(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_IoReader::reg_IoReader_void() {
  int epidx = CkRegisterEp("IoReader()",
      reinterpret_cast<CkCallFnPtr>(_call_IoReader_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoReader::_call_IoReader_void(void* impl_msg, void* impl_obj_void)
{
  IoReader* impl_obj = static_cast<IoReader*>(impl_obj_void);
  new (impl_obj_void) IoReader();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoReader(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_IoReader::reg_IoReader_CkMigrateMessage() {
  int epidx = CkRegisterEp("IoReader(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_IoReader_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoReader::_call_IoReader_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  call_migration_constructor<IoReader> c = impl_obj_void;
  c((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoReader();
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoReader(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_IoReader::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterArrayDimensions(__idx, 1);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
  // REG: IoReader();
  idx_IoReader_void();
  CkRegisterDefaultCtor(__idx, idx_IoReader_void());

  // REG: IoReader(CkMigrateMessage* impl_msg);
  idx_IoReader_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_IoReader_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: array IoWriter: ArrayElement{
IoWriter();
IoWriter(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_IoWriter::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CProxySection_IoWriter::contribute(CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, userData, fragSize);
}

void CProxySection_IoWriter::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, userData, fragSize);
}

template <typename T>
void CProxySection_IoWriter::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, userData, fragSize);
}

void CProxySection_IoWriter::contribute(CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, cb, userData, fragSize);
}

void CProxySection_IoWriter::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, cb, userData, fragSize);
}

template <typename T>
void CProxySection_IoWriter::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, cb, userData, fragSize);
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoWriter();
 */
void CProxyElement_IoWriter::insert(int onPE, const CkEntryOptions *impl_e_opts)
{ 
   void *impl_msg = CkAllocSysMsg(impl_e_opts);
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_IoWriter::idx_IoWriter_void(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoWriter(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoWriter();
 */
CkArrayID CProxy_IoWriter::ckNew(const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoWriter::idx_IoWriter_void(), opts);
  return gId;
}
void CProxy_IoWriter::ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoWriter::idx_IoWriter_void(), _ck_array_creation_cb, opts, impl_msg);
}
CkArrayID CProxy_IoWriter::ckNew(const int s1, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoWriter::idx_IoWriter_void(), opts);
  return gId;
}
void CProxy_IoWriter::ckNew(const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoWriter::idx_IoWriter_void(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_IoWriter::reg_IoWriter_void() {
  int epidx = CkRegisterEp("IoWriter()",
      reinterpret_cast<CkCallFnPtr>(_call_IoWriter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoWriter::_call_IoWriter_void(void* impl_msg, void* impl_obj_void)
{
  IoWriter* impl_obj = static_cast<IoWriter*>(impl_obj_void);
  new (impl_obj_void) IoWriter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoWriter(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_IoWriter::reg_IoWriter_CkMigrateMessage() {
  int epidx = CkRegisterEp("IoWriter(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_IoWriter_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoWriter::_call_IoWriter_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  call_migration_constructor<IoWriter> c = impl_obj_void;
  c((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoWriter();
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoWriter(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_IoWriter::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterArrayDimensions(__idx, 1);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
  // REG: IoWriter();
  idx_IoWriter_void();
  CkRegisterDefaultCtor(__idx, idx_IoWriter_void());

  // REG: IoWriter(CkMigrateMessage* impl_msg);
  idx_IoWriter_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_IoWriter_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
void _registermesh(void)
{
  static int _done = 0; if(_done) return; _done = 1;
  _registerInitCall(register_reduce_performance,1);

  _registerInitCall(register_reduce_method_debug,1);

  _registerInitCall(register_sum_long_double,1);

  _registerInitCall(register_sum_long_double_2,1);

  _registerInitCall(register_sum_long_double_3,1);

  _registerInitCall(register_sum_long_double_n,1);

  _registerInitCall(mutex_init_hierarchy,1);

  _registerInitCall(mutex_init_initial_value,1);

  _registerInitCall(mutex_init_field_face,1);

  CkRegisterReadonly("MsgCoarsen::counter","int",sizeof(MsgCoarsen::counter),(void *) &MsgCoarsen::counter,__xlater_roPup_MsgCoarsen_QColon__QColon_counter);

  CkRegisterReadonly("MsgAdapt::counter","int",sizeof(MsgAdapt::counter),(void *) &MsgAdapt::counter,__xlater_roPup_MsgAdapt_QColon__QColon_counter);

  CkRegisterReadonly("MsgInitial::counter","int",sizeof(MsgInitial::counter),(void *) &MsgInitial::counter,__xlater_roPup_MsgInitial_QColon__QColon_counter);

  CkRegisterReadonly("MsgOutput::counter","int",sizeof(MsgOutput::counter),(void *) &MsgOutput::counter,__xlater_roPup_MsgOutput_QColon__QColon_counter);

  CkRegisterReadonly("MsgRefine::counter","int",sizeof(MsgRefine::counter),(void *) &MsgRefine::counter,__xlater_roPup_MsgRefine_QColon__QColon_counter);

  CkRegisterReadonly("MsgRefresh::counter","int",sizeof(MsgRefresh::counter),(void *) &MsgRefresh::counter,__xlater_roPup_MsgRefresh_QColon__QColon_counter);

  CkRegisterReadonly("DataMsg::counter","int",sizeof(DataMsg::counter),(void *) &DataMsg::counter,__xlater_roPup_DataMsg_QColon__QColon_counter);

  CkRegisterReadonly("FieldFace::counter","int",sizeof(FieldFace::counter),(void *) &FieldFace::counter,__xlater_roPup_FieldFace_QColon__QColon_counter);

  CkRegisterReadonly("ParticleData::counter","int",sizeof(ParticleData::counter),(void *) &ParticleData::counter,__xlater_roPup_ParticleData_QColon__QColon_counter);

  CkRegisterReadonly("InitialTrace::id0_","int",sizeof(InitialTrace::id0_),(void *) &InitialTrace::id0_,__xlater_roPup_InitialTrace_QColon__QColon_id0_);

  CkRegisterReadonly("Method::courant_global","double",sizeof(Method::courant_global),(void *) &Method::courant_global,__xlater_roPup_Method_QColon__QColon_courant_global);

  CkRegisterReadonly("g_config","Config",sizeof(g_config),(void *) &g_config,__xlater_roPup_g_config);

  CkRegisterReadonly("g_parameters","Parameters",sizeof(g_parameters),(void *) &g_parameters,__xlater_roPup_g_parameters);

      PUPable_reg(BoundaryPeriodic);

      PUPable_reg(BoundaryValue);

      PUPable_reg(ColormapRGB);

      PUPable_reg(Config);

      PUPable_reg(Factory);

      PUPable_reg(Initial);

      PUPable_reg(InitialTrace);

      PUPable_reg(InitialValue);

      PUPable_reg(Io);

      PUPable_reg(IoBlock);

      PUPable_reg(IoFieldData);

      PUPable_reg(IoHierarchy);

      PUPable_reg(IoParticleData);

      PUPable_reg(IoSimulation);

      PUPable_reg(ItIndexList);

      PUPable_reg(ItIndexRange);

      PUPable_reg(MaskExpr);

      PUPable_reg(MaskPng);

      PUPable_reg(MethodCloseFiles);

      PUPable_reg(MethodDebug);

      PUPable_reg(MethodFluxCorrect);

      PUPable_reg(MethodNull);

      PUPable_reg(MethodOrderMorton);

      PUPable_reg(MethodOutput);

      PUPable_reg(MethodRefresh);

      PUPable_reg(MethodTrace);

      PUPable_reg(ObjectSphere);

      PUPable_reg(OutputCheckpoint);

      PUPable_reg(OutputData);

      PUPable_reg(OutputImage);

      PUPable_reg(Physics);

      PUPable_reg(Problem);

      PUPable_reg(ProlongInject);

      PUPable_reg(ProlongLinear);

      PUPable_reg(Refine);

      PUPable_reg(RefineDensity);

      PUPable_reg(RefineMask);

      PUPable_reg(RefineParticleCount);

      PUPable_reg(RefineShear);

      PUPable_reg(RefineSlope);

      PUPable_reg(Refresh);

      PUPable_reg(RestrictLinear);

      PUPable_reg(ScheduleInterval);

      PUPable_reg(ScheduleList);

      PUPable_reg(SolverNull);

      PUPable_reg(Stopping);

      PUPable_reg(Units);

/* REG: message FieldMsg{
char a[];
}
;
*/
CMessage_FieldMsg::__register("FieldMsg", sizeof(FieldMsg),(CkPackFnPtr) FieldMsg::pack,(CkUnpackFnPtr) FieldMsg::unpack);

/* REG: message MsgAdapt;
*/
CMessage_MsgAdapt::__register("MsgAdapt", sizeof(MsgAdapt),(CkPackFnPtr) MsgAdapt::pack,(CkUnpackFnPtr) MsgAdapt::unpack);

/* REG: message MsgCoarsen;
*/
CMessage_MsgCoarsen::__register("MsgCoarsen", sizeof(MsgCoarsen),(CkPackFnPtr) MsgCoarsen::pack,(CkUnpackFnPtr) MsgCoarsen::unpack);

/* REG: message MsgInitial;
*/
CMessage_MsgInitial::__register("MsgInitial", sizeof(MsgInitial),(CkPackFnPtr) MsgInitial::pack,(CkUnpackFnPtr) MsgInitial::unpack);

/* REG: message MsgOutput;
*/
CMessage_MsgOutput::__register("MsgOutput", sizeof(MsgOutput),(CkPackFnPtr) MsgOutput::pack,(CkUnpackFnPtr) MsgOutput::unpack);

/* REG: message MsgRefine;
*/
CMessage_MsgRefine::__register("MsgRefine", sizeof(MsgRefine),(CkPackFnPtr) MsgRefine::pack,(CkUnpackFnPtr) MsgRefine::unpack);

/* REG: message MsgRefresh;
*/
CMessage_MsgRefresh::__register("MsgRefresh", sizeof(MsgRefresh),(CkPackFnPtr) MsgRefresh::pack,(CkUnpackFnPtr) MsgRefresh::unpack);

/* REG: array Block: ArrayElement{
Block(const process_type &ip_source, const MsgType &msg_type);
void p_set_msg_refine(MsgRefine* impl_msg);
Block();
void p_initial_exit();
void r_end_initialize(CkReductionMsg* impl_msg);
void r_initial_new_next(CkReductionMsg* impl_msg);
void r_initial_new_continue(CkReductionMsg* impl_msg);
void p_output_enter();
void r_output_enter(CkReductionMsg* impl_msg);
void p_output_end();
void p_output_exit();
void r_output_exit(CkReductionMsg* impl_msg);
void p_output_write(int index_output, int step);
void p_compute_enter();
void p_compute_continue();
void r_compute_continue(CkReductionMsg* impl_msg);
void p_compute_exit();
void r_compute_exit(CkReductionMsg* impl_msg);
void p_method_flux_correct_refresh();
void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
void r_method_order_morton_continue(CkReductionMsg* impl_msg);
void r_method_order_morton_complete(CkReductionMsg* impl_msg);
void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
void p_method_order_morton_index(int index, int count);
void p_method_output_next(MsgOutput* impl_msg);
void p_method_output_write(MsgOutput* impl_msg);
void r_method_output_continue(CkReductionMsg* impl_msg);
void r_method_output_done(CkReductionMsg* impl_msg);
void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
void p_stopping_enter();
void r_stopping_enter(CkReductionMsg* impl_msg);
void p_stopping_load_balance();
void r_stopping_load_balance(CkReductionMsg* impl_msg);
void p_stopping_exit();
void r_stopping_exit(CkReductionMsg* impl_msg);
void p_exit();
void r_exit(CkReductionMsg* impl_msg);
void p_control_sync_count(int entry_point, int id, int count);
void p_adapt_enter();
void r_adapt_enter(CkReductionMsg* impl_msg);
void p_adapt_end();
void p_adapt_update();
void r_adapt_next(CkReductionMsg* impl_msg);
void p_adapt_called();
void p_adapt_exit();
void p_adapt_delete();
void p_adapt_recv_level(MsgAdapt* impl_msg);
void p_adapt_recv_child(MsgCoarsen* impl_msg);
void r_restart_enter(CkReductionMsg* impl_msg);
void p_refresh_recv(MsgRefresh* impl_msg);
void p_refresh_child(int n, const char *a, const int *ic3);
Block(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_Block::__register("Block", sizeof(Block));

/* REG: array IoReader: ArrayElement{
IoReader();
IoReader(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_IoReader::__register("IoReader", sizeof(IoReader));

/* REG: array IoWriter: ArrayElement{
IoWriter();
IoWriter(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_IoWriter::__register("IoWriter", sizeof(IoWriter));

}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_Block::virtual_pup(PUP::er &p) {
    recursive_pup<Block>(dynamic_cast<Block*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_IoReader::virtual_pup(PUP::er &p) {
    recursive_pup<IoReader>(dynamic_cast<IoReader*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_IoWriter::virtual_pup(PUP::er &p) {
    recursive_pup<IoWriter>(dynamic_cast<IoWriter*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
