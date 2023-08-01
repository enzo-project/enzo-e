



















































































































/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_get_msg_refine_3_closure : public SDAG::Closure {
            Index index;


      p_get_msg_refine_3_closure() {
        init();
      }
      p_get_msg_refine_3_closure(CkMigrateMessage*) {
        init();
      }
            Index & getP0() { return index;}
      void pup(PUP::er& __p) {
        __p | index;
        packClosure(__p);
      }
      virtual ~p_get_msg_refine_3_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_get_msg_refine_3_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_get_msg_check_4_closure : public SDAG::Closure {
            Index index;


      p_get_msg_check_4_closure() {
        init();
      }
      p_get_msg_check_4_closure(CkMigrateMessage*) {
        init();
      }
            Index & getP0() { return index;}
      void pup(PUP::er& __p) {
        __p | index;
        packClosure(__p);
      }
      virtual ~p_get_msg_check_4_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_get_msg_check_4_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_method_balance_check_6_closure : public SDAG::Closure {
      

      p_method_balance_check_6_closure() {
        init();
      }
      p_method_balance_check_6_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_balance_check_6_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_balance_check_6_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_check_done_8_closure : public SDAG::Closure {
      

      p_check_done_8_closure() {
        init();
      }
      p_check_done_8_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_check_done_8_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_check_done_8_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_set_io_writer_9_closure : public SDAG::Closure {
            CProxy_IoEnzoWriter proxy;


      p_set_io_writer_9_closure() {
        init();
      }
      p_set_io_writer_9_closure(CkMigrateMessage*) {
        init();
      }
            CProxy_IoEnzoWriter & getP0() { return proxy;}
      void pup(PUP::er& __p) {
        __p | proxy;
        packClosure(__p);
      }
      virtual ~p_set_io_writer_9_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_set_io_writer_9_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_infer_set_array_count_10_closure : public SDAG::Closure {
            int count;


      p_infer_set_array_count_10_closure() {
        init();
      }
      p_infer_set_array_count_10_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return count;}
      void pup(PUP::er& __p) {
        __p | count;
        packClosure(__p);
      }
      virtual ~p_infer_set_array_count_10_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_infer_set_array_count_10_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_infer_array_created_11_closure : public SDAG::Closure {
      

      p_infer_array_created_11_closure() {
        init();
      }
      p_infer_array_created_11_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_infer_array_created_11_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_infer_array_created_11_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_infer_done_12_closure : public SDAG::Closure {
      

      p_infer_done_12_closure() {
        init();
      }
      p_infer_done_12_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_infer_done_12_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_infer_done_12_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_set_io_reader_13_closure : public SDAG::Closure {
            CProxy_IoEnzoReader proxy;


      p_set_io_reader_13_closure() {
        init();
      }
      p_set_io_reader_13_closure(CkMigrateMessage*) {
        init();
      }
            CProxy_IoEnzoReader & getP0() { return proxy;}
      void pup(PUP::er& __p) {
        __p | proxy;
        packClosure(__p);
      }
      virtual ~p_set_io_reader_13_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_set_io_reader_13_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_io_reader_created_14_closure : public SDAG::Closure {
      

      p_io_reader_created_14_closure() {
        init();
      }
      p_io_reader_created_14_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_io_reader_created_14_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_io_reader_created_14_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_restart_next_level_15_closure : public SDAG::Closure {
      

      p_restart_next_level_15_closure() {
        init();
      }
      p_restart_next_level_15_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_restart_next_level_15_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_restart_next_level_15_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_restart_level_created_16_closure : public SDAG::Closure {
      

      p_restart_level_created_16_closure() {
        init();
      }
      p_restart_level_created_16_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_restart_level_created_16_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_restart_level_created_16_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_set_level_array_17_closure : public SDAG::Closure {
            CProxy_EnzoLevelArray proxy;


      p_set_level_array_17_closure() {
        init();
      }
      p_set_level_array_17_closure(CkMigrateMessage*) {
        init();
      }
            CProxy_EnzoLevelArray & getP0() { return proxy;}
      void pup(PUP::er& __p) {
        __p | proxy;
        packClosure(__p);
      }
      virtual ~p_set_level_array_17_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_set_level_array_17_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_fbnet_concatenate_sphere_lists_18_closure : public SDAG::Closure {
      

      p_fbnet_concatenate_sphere_lists_18_closure() {
        init();
      }
      p_fbnet_concatenate_sphere_lists_18_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_fbnet_concatenate_sphere_lists_18_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_fbnet_concatenate_sphere_lists_18_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoSimulation::p_fbnet_done_19_closure : public SDAG::Closure {
      

      p_fbnet_done_19_closure() {
        init();
      }
      p_fbnet_done_19_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_fbnet_done_19_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_fbnet_done_19_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_feedback_starss_end_5_closure : public SDAG::Closure {
      

      p_method_feedback_starss_end_5_closure() {
        init();
      }
      p_method_feedback_starss_end_5_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_feedback_starss_end_5_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_feedback_starss_end_5_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_m1_closure_solve_transport_eqn_6_closure : public SDAG::Closure {
      

      p_method_m1_closure_solve_transport_eqn_6_closure() {
        init();
      }
      p_method_m1_closure_solve_transport_eqn_6_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_m1_closure_solve_transport_eqn_6_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_m1_closure_solve_transport_eqn_6_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_balance_migrate_10_closure : public SDAG::Closure {
      

      p_method_balance_migrate_10_closure() {
        init();
      }
      p_method_balance_migrate_10_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_balance_migrate_10_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_balance_migrate_10_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_balance_done_11_closure : public SDAG::Closure {
      

      p_method_balance_done_11_closure() {
        init();
      }
      p_method_balance_done_11_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_balance_done_11_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_balance_done_11_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_gravity_continue_12_closure : public SDAG::Closure {
      

      p_method_gravity_continue_12_closure() {
        init();
      }
      p_method_gravity_continue_12_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_gravity_continue_12_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_gravity_continue_12_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_gravity_end_13_closure : public SDAG::Closure {
      

      p_method_gravity_end_13_closure() {
        init();
      }
      p_method_gravity_end_13_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_gravity_end_13_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_gravity_end_13_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_infer_merge_masks_14_closure : public SDAG::Closure {
            int n;
            char *mask;
            int *ic3;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_method_infer_merge_masks_14_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_method_infer_merge_masks_14_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return n;}
            char *& getP1() { return mask;}
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
  int impl_off_mask, impl_cnt_mask;
  implP|impl_off_mask;
  implP|impl_cnt_mask;
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
          impl_buf+=CK_ALIGN(implP.size(),16);
          mask = (char *)(impl_buf+impl_off_mask);
          ic3 = (int *)(impl_buf+impl_off_ic3);
        }
      }
      virtual ~p_method_infer_merge_masks_14_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_method_infer_merge_masks_14_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_infer_count_arrays_15_closure : public SDAG::Closure {
            int count;


      p_method_infer_count_arrays_15_closure() {
        init();
      }
      p_method_infer_count_arrays_15_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return count;}
      void pup(PUP::er& __p) {
        __p | count;
        packClosure(__p);
      }
      virtual ~p_method_infer_count_arrays_15_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_infer_count_arrays_15_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_infer_request_data_16_closure : public SDAG::Closure {
            int *il3;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_method_infer_request_data_16_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_method_infer_request_data_16_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int *& getP0() { return il3;}
      void pup(PUP::er& __p) {
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  int impl_off_il3, impl_cnt_il3;
  implP|impl_off_il3;
  implP|impl_cnt_il3;
          impl_buf+=CK_ALIGN(implP.size(),16);
          il3 = (int *)(impl_buf+impl_off_il3);
        }
      }
      virtual ~p_method_infer_request_data_16_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_method_infer_request_data_16_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_infer_update_17_closure : public SDAG::Closure {
            int n;
            char *buffer;
            int *il3;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_method_infer_update_17_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_method_infer_update_17_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return n;}
            char *& getP1() { return buffer;}
            int *& getP2() { return il3;}
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
  int impl_off_buffer, impl_cnt_buffer;
  implP|impl_off_buffer;
  implP|impl_cnt_buffer;
  int impl_off_il3, impl_cnt_il3;
  implP|impl_off_il3;
  implP|impl_cnt_il3;
          impl_buf+=CK_ALIGN(implP.size(),16);
          buffer = (char *)(impl_buf+impl_off_buffer);
          il3 = (int *)(impl_buf+impl_off_il3);
        }
      }
      virtual ~p_method_infer_update_17_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_method_infer_update_17_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_infer_exit_18_closure : public SDAG::Closure {
      

      p_method_infer_exit_18_closure() {
        init();
      }
      p_method_infer_exit_18_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_infer_exit_18_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_infer_exit_18_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_fbnet_exit_20_closure : public SDAG::Closure {
      

      p_method_fbnet_exit_20_closure() {
        init();
      }
      p_method_fbnet_exit_20_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_fbnet_exit_20_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_fbnet_exit_20_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_check_write_first_21_closure : public SDAG::Closure {
            int num_files;
            std::string ordering;
            std::string name_dir;


      p_check_write_first_21_closure() {
        init();
      }
      p_check_write_first_21_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return num_files;}
            std::string & getP1() { return ordering;}
            std::string & getP2() { return name_dir;}
      void pup(PUP::er& __p) {
        __p | num_files;
        __p | ordering;
        __p | name_dir;
        packClosure(__p);
      }
      virtual ~p_check_write_first_21_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_check_write_first_21_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_check_write_next_22_closure : public SDAG::Closure {
            int num_files;
            std::string ordering;


      p_check_write_next_22_closure() {
        init();
      }
      p_check_write_next_22_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return num_files;}
            std::string & getP1() { return ordering;}
      void pup(PUP::er& __p) {
        __p | num_files;
        __p | ordering;
        packClosure(__p);
      }
      virtual ~p_check_write_next_22_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_check_write_next_22_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_check_done_23_closure : public SDAG::Closure {
      

      p_check_done_23_closure() {
        init();
      }
      p_check_done_23_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_check_done_23_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_check_done_23_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_restart_refine_24_closure : public SDAG::Closure {
            int *ic3;
            int io_reader;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_restart_refine_24_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_restart_refine_24_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int *& getP0() { return ic3;}
            int & getP1() { return io_reader;}
      void pup(PUP::er& __p) {
        __p | io_reader;
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
  PUP::detail::TemporaryObjectHolder<int> io_reader;
  implP|io_reader;
          impl_buf+=CK_ALIGN(implP.size(),16);
          ic3 = (int *)(impl_buf+impl_off_ic3);
        }
      }
      virtual ~p_restart_refine_24_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_restart_refine_24_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_restart_done_26_closure : public SDAG::Closure {
      

      p_restart_done_26_closure() {
        init();
      }
      p_restart_done_26_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_restart_done_26_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_restart_done_26_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_method_accretion_end_27_closure : public SDAG::Closure {
      

      p_method_accretion_end_27_closure() {
        init();
      }
      p_method_accretion_end_27_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_method_accretion_end_27_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_method_accretion_end_27_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_cg_matvec_28_closure : public SDAG::Closure {
      

      p_solver_cg_matvec_28_closure() {
        init();
      }
      p_solver_cg_matvec_28_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_cg_matvec_28_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_cg_matvec_28_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_cg_loop_2_32_closure : public SDAG::Closure {
      

      p_solver_cg_loop_2_32_closure() {
        init();
      }
      p_solver_cg_loop_2_32_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_cg_loop_2_32_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_cg_loop_2_32_closure));
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
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_bicgstab_loop_2_41_closure : public SDAG::Closure {
      

      p_solver_bicgstab_loop_2_41_closure() {
        init();
      }
      p_solver_bicgstab_loop_2_41_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_bicgstab_loop_2_41_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_bicgstab_loop_2_41_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_bicgstab_loop_3_42_closure : public SDAG::Closure {
      

      p_solver_bicgstab_loop_3_42_closure() {
        init();
      }
      p_solver_bicgstab_loop_3_42_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_bicgstab_loop_3_42_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_bicgstab_loop_3_42_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_bicgstab_loop_8_43_closure : public SDAG::Closure {
      

      p_solver_bicgstab_loop_8_43_closure() {
        init();
      }
      p_solver_bicgstab_loop_8_43_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_bicgstab_loop_8_43_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_bicgstab_loop_8_43_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_bicgstab_loop_9_44_closure : public SDAG::Closure {
      

      p_solver_bicgstab_loop_9_44_closure() {
        init();
      }
      p_solver_bicgstab_loop_9_44_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_bicgstab_loop_9_44_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_bicgstab_loop_9_44_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_dot_recv_parent_45_closure : public SDAG::Closure {
            int n;
            long double *dot;
            std::vector<int> isa;
            int i_function;
            int iter;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_dot_recv_parent_45_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_dot_recv_parent_45_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return n;}
            long double *& getP1() { return dot;}
            std::vector<int> & getP2() { return isa;}
            int & getP3() { return i_function;}
            int & getP4() { return iter;}
      void pup(PUP::er& __p) {
        __p | n;
        __p | isa;
        __p | i_function;
        __p | iter;
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
  int impl_off_dot, impl_cnt_dot;
  implP|impl_off_dot;
  implP|impl_cnt_dot;
  PUP::detail::TemporaryObjectHolder<std::vector<int>> isa;
  implP|isa;
  PUP::detail::TemporaryObjectHolder<int> i_function;
  implP|i_function;
  PUP::detail::TemporaryObjectHolder<int> iter;
  implP|iter;
          impl_buf+=CK_ALIGN(implP.size(),16);
          dot = (long double *)(impl_buf+impl_off_dot);
        }
      }
      virtual ~p_dot_recv_parent_45_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_dot_recv_parent_45_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_dot_recv_children_46_closure : public SDAG::Closure {
            int n;
            long double *dot;
            std::vector<int> isa;
            int i_function;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_dot_recv_children_46_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_dot_recv_children_46_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return n;}
            long double *& getP1() { return dot;}
            std::vector<int> & getP2() { return isa;}
            int & getP3() { return i_function;}
      void pup(PUP::er& __p) {
        __p | n;
        __p | isa;
        __p | i_function;
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
  int impl_off_dot, impl_cnt_dot;
  implP|impl_off_dot;
  implP|impl_cnt_dot;
  PUP::detail::TemporaryObjectHolder<std::vector<int>> isa;
  implP|isa;
  PUP::detail::TemporaryObjectHolder<int> i_function;
  implP|i_function;
          impl_buf+=CK_ALIGN(implP.size(),16);
          dot = (long double *)(impl_buf+impl_off_dot);
        }
      }
      virtual ~p_dot_recv_children_46_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_dot_recv_children_46_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_dd_solve_coarse_49_closure : public SDAG::Closure {
      

      p_solver_dd_solve_coarse_49_closure() {
        init();
      }
      p_solver_dd_solve_coarse_49_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_dd_solve_coarse_49_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_dd_solve_coarse_49_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_dd_solve_domain_50_closure : public SDAG::Closure {
      

      p_solver_dd_solve_domain_50_closure() {
        init();
      }
      p_solver_dd_solve_domain_50_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_dd_solve_domain_50_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_dd_solve_domain_50_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_dd_last_smooth_51_closure : public SDAG::Closure {
      

      p_solver_dd_last_smooth_51_closure() {
        init();
      }
      p_solver_dd_last_smooth_51_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_dd_last_smooth_51_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_dd_last_smooth_51_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_jacobi_continue_54_closure : public SDAG::Closure {
      

      p_solver_jacobi_continue_54_closure() {
        init();
      }
      p_solver_jacobi_continue_54_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_jacobi_continue_54_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_jacobi_continue_54_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_mg0_restrict_55_closure : public SDAG::Closure {
      

      p_solver_mg0_restrict_55_closure() {
        init();
      }
      p_solver_mg0_restrict_55_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_mg0_restrict_55_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_mg0_restrict_55_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_mg0_solve_coarse_56_closure : public SDAG::Closure {
      

      p_solver_mg0_solve_coarse_56_closure() {
        init();
      }
      p_solver_mg0_solve_coarse_56_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_mg0_solve_coarse_56_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_mg0_solve_coarse_56_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_mg0_post_smooth_57_closure : public SDAG::Closure {
      

      p_solver_mg0_post_smooth_57_closure() {
        init();
      }
      p_solver_mg0_post_smooth_57_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_mg0_post_smooth_57_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_mg0_post_smooth_57_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoBlock::p_solver_mg0_last_smooth_58_closure : public SDAG::Closure {
      

      p_solver_mg0_last_smooth_58_closure() {
        init();
      }
      p_solver_mg0_last_smooth_58_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_solver_mg0_last_smooth_58_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_solver_mg0_last_smooth_58_closure));
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


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_IoEnzoReader::p_init_root_2_closure : public SDAG::Closure {
            std::string impl_noname_0;
            std::string impl_noname_1;
            int level;


      p_init_root_2_closure() {
        init();
      }
      p_init_root_2_closure(CkMigrateMessage*) {
        init();
      }
            std::string & getP0() { return impl_noname_0;}
            std::string & getP1() { return impl_noname_1;}
            int & getP2() { return level;}
      void pup(PUP::er& __p) {
        __p | impl_noname_0;
        __p | impl_noname_1;
        __p | level;
        packClosure(__p);
      }
      virtual ~p_init_root_2_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_init_root_2_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_IoEnzoReader::p_create_level_3_closure : public SDAG::Closure {
            int level;


      p_create_level_3_closure() {
        init();
      }
      p_create_level_3_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return level;}
      void pup(PUP::er& __p) {
        __p | level;
        packClosure(__p);
      }
      virtual ~p_create_level_3_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_create_level_3_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_IoEnzoReader::p_init_level_4_closure : public SDAG::Closure {
            int level;


      p_init_level_4_closure() {
        init();
      }
      p_init_level_4_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return level;}
      void pup(PUP::er& __p) {
        __p | level;
        packClosure(__p);
      }
      virtual ~p_init_level_4_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_init_level_4_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_IoEnzoReader::p_block_created_5_closure : public SDAG::Closure {
      

      p_block_created_5_closure() {
        init();
      }
      p_block_created_5_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_block_created_5_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_block_created_5_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_IoEnzoReader::p_block_ready_6_closure : public SDAG::Closure {
      

      p_block_ready_6_closure() {
        init();
      }
      p_block_ready_6_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_block_ready_6_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_block_ready_6_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoLevelArray::p_request_data_2_closure : public SDAG::Closure {
      

      p_request_data_2_closure() {
        init();
      }
      p_request_data_2_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_request_data_2_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_request_data_2_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoLevelArray::p_transfer_data_3_closure : public SDAG::Closure {
            Index impl_noname_3;
            int nf;
            enzo_float *field_data;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_transfer_data_3_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_transfer_data_3_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            Index & getP0() { return impl_noname_3;}
            int & getP1() { return nf;}
            enzo_float *& getP2() { return field_data;}
      void pup(PUP::er& __p) {
        __p | impl_noname_3;
        __p | nf;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_3;
  implP|impl_noname_3;
  PUP::detail::TemporaryObjectHolder<int> nf;
  implP|nf;
  int impl_off_field_data, impl_cnt_field_data;
  implP|impl_off_field_data;
  implP|impl_cnt_field_data;
          impl_buf+=CK_ALIGN(implP.size(),16);
          field_data = (enzo_float *)(impl_buf+impl_off_field_data);
        }
      }
      virtual ~p_transfer_data_3_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_transfer_data_3_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_EnzoLevelArray::p_done_4_closure : public SDAG::Closure {
            Index impl_noname_4;


      p_done_4_closure() {
        init();
      }
      p_done_4_closure(CkMigrateMessage*) {
        init();
      }
            Index & getP0() { return impl_noname_4;}
      void pup(PUP::er& __p) {
        __p | impl_noname_4;
        packClosure(__p);
      }
      virtual ~p_done_4_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_done_4_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */





/* DEFS: readonly int EnzoMsgCheck::counter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoMsgCheck_QColon__QColon_counter(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoMsgCheck::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoMsgCheck::counter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoMsgCheck::counter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly EnzoConfig g_enzo_config;
 */
extern EnzoConfig g_enzo_config;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_g_enzo_config(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|g_enzo_config;
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::UseMinimumPressureSupport[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_UseMinimumPressureSupport(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoBlock::UseMinimumPressureSupport[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::UseMinimumPressureSupport[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::UseMinimumPressureSupport[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly enzo_float EnzoBlock::MinimumPressureSupportParameter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_MinimumPressureSupportParameter(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(enzo_float) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& EnzoBlock::MinimumPressureSupportParameter[0], ((CONFIG_NODE_SIZE) * sizeof(enzo_float)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::MinimumPressureSupportParameter[0], ((CONFIG_NODE_SIZE) * sizeof(enzo_float)), regMode);
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
  _impl_p(&EnzoBlock::MinimumPressureSupportParameter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::MultiSpecies[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_MultiSpecies(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoBlock::MultiSpecies[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::MultiSpecies[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::MultiSpecies[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::PressureFree[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_PressureFree(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoBlock::PressureFree[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::PressureFree[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::PressureFree[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly enzo_float EnzoBlock::GravitationalConstant[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_GravitationalConstant(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(enzo_float) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& EnzoBlock::GravitationalConstant[0], ((CONFIG_NODE_SIZE) * sizeof(enzo_float)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::GravitationalConstant[0], ((CONFIG_NODE_SIZE) * sizeof(enzo_float)), regMode);
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
  _impl_p(&EnzoBlock::GravitationalConstant[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::ProblemType[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_ProblemType(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoBlock::ProblemType[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::ProblemType[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::ProblemType[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::PPMFlatteningParameter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_PPMFlatteningParameter(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoBlock::PPMFlatteningParameter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::PPMFlatteningParameter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::PPMFlatteningParameter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::PPMDiffusionParameter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_PPMDiffusionParameter(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoBlock::PPMDiffusionParameter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::PPMDiffusionParameter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::PPMDiffusionParameter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::PPMSteepeningParameter[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_PPMSteepeningParameter(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoBlock::PPMSteepeningParameter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::PPMSteepeningParameter[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::PPMSteepeningParameter[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly enzo_float EnzoBlock::InitialRedshift[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_InitialRedshift(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(enzo_float) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& EnzoBlock::InitialRedshift[0], ((CONFIG_NODE_SIZE) * sizeof(enzo_float)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::InitialRedshift[0], ((CONFIG_NODE_SIZE) * sizeof(enzo_float)), regMode);
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
  _impl_p(&EnzoBlock::InitialRedshift[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly enzo_float EnzoBlock::InitialTimeInCodeUnits[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_InitialTimeInCodeUnits(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE) * sizeof(enzo_float) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& EnzoBlock::InitialTimeInCodeUnits[0], ((CONFIG_NODE_SIZE) * sizeof(enzo_float)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::InitialTimeInCodeUnits[0], ((CONFIG_NODE_SIZE) * sizeof(enzo_float)), regMode);
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
  _impl_p(&EnzoBlock::InitialTimeInCodeUnits[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly enzo_float EnzoBlock::DomainLeftEdge[CONFIG_NODE_SIZE_3];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_DomainLeftEdge(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE_3) * sizeof(enzo_float) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& EnzoBlock::DomainLeftEdge[0], ((CONFIG_NODE_SIZE_3) * sizeof(enzo_float)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::DomainLeftEdge[0], ((CONFIG_NODE_SIZE_3) * sizeof(enzo_float)), regMode);
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
  _impl_p(&EnzoBlock::DomainLeftEdge[0], ((CONFIG_NODE_SIZE_3)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly enzo_float EnzoBlock::DomainRightEdge[CONFIG_NODE_SIZE_3];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_DomainRightEdge(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE_3) * sizeof(enzo_float) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& EnzoBlock::DomainRightEdge[0], ((CONFIG_NODE_SIZE_3) * sizeof(enzo_float)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::DomainRightEdge[0], ((CONFIG_NODE_SIZE_3) * sizeof(enzo_float)), regMode);
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
  _impl_p(&EnzoBlock::DomainRightEdge[0], ((CONFIG_NODE_SIZE_3)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::GridRank[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_GridRank(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoBlock::GridRank[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::GridRank[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::GridRank[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::ghost_depth[CONFIG_NODE_SIZE_3];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_ghost_depth(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
#if !CMK_ONESIDED_RO_DISABLE && CMK_ONESIDED_IMPL
  if((CONFIG_NODE_SIZE_3) * sizeof(int) >= CMK_ONESIDED_RO_THRESHOLD) {
    if(_impl_p.isSizing()) {
      CkNcpyBuffer myBuffer;
      _impl_p|myBuffer;
      readonlyUpdateNumops();
    }
    if(_impl_p.isPacking()) {
    int regMode = CK_BUFFER_REG;
    if(CkNumNodes() == 1)
      regMode = CK_BUFFER_UNREG;
    CkNcpyBuffer myBuffer(& EnzoBlock::ghost_depth[0], ((CONFIG_NODE_SIZE_3) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::ghost_depth[0], ((CONFIG_NODE_SIZE_3) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::ghost_depth[0], ((CONFIG_NODE_SIZE_3)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly int EnzoBlock::NumberOfBaryonFields[CONFIG_NODE_SIZE];
 */
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_EnzoBlock_QColon__QColon_NumberOfBaryonFields(void *_impl_pup_er) {
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
    CkNcpyBuffer myBuffer(& EnzoBlock::NumberOfBaryonFields[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
    CkNcpyBuffer myBuffer(& EnzoBlock::NumberOfBaryonFields[0], ((CONFIG_NODE_SIZE) * sizeof(int)), regMode);
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
  _impl_p(&EnzoBlock::NumberOfBaryonFields[0], ((CONFIG_NODE_SIZE)) );
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoBoundary)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoConfig)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoFactory)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialGrackleTest)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialAccretionTest)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialBBTest)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialBCenter)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialBurkertBodenheimer)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialCloud)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialCollapse)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialCosmology)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialFeedbackTest)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialHdf5)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialImplosion2)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialInclinedWave)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialIsolatedGalaxy)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialMergeSinksTest)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialMusic)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialPm)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialPpmlTest)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialSedovArray2)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialSedovArray3)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialSedovRandom)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialShockTube)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialShuCollapse)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialSoup)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoInitialTurbulence)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(IoEnzoBlock)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoRefineMass)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoRefineParticleMass)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoRefineShock)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoComputeCoolingTime)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoComputePressure)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoComputeTemperature)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoComputeAcceleration)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoComputeCicInterp)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMatrixLaplace)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMatrixDiagonal)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMatrixIdentity)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodGrackle)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodAccretion)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodBackgroundAcceleration)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodBondiHoyleAccretion)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodCheck)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodComovingExpansion)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodCosmology)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodDistributedFeedback)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodFeedback)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodFeedbackSTARSS)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodM1Closure)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodFluxAccretion)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodGravity)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodHeat)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodInference)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodFBNetDeposit)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodMergeSinks)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodMHDVlct)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodPmDeposit)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodPmUpdate)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodPpm)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodPpml)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodBalance)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodSinkMaker)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodThresholdAccretion)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodTurbulence)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoIntegrationQuanUpdate)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoReconstructorNN)
#endif /* CK_TEMPLATES_ONLY */

#ifdef CK_TEMPLATES_ONLY
  #define _CHARMXI_CLASS_NAME EnzoReconstructorPLM<PLM_EnzoRKLimiter>
  PUPable_def_template(_CHARMXI_CLASS_NAME)
  #undef _CHARMXI_CLASS_NAME
#endif /* CK_TEMPLATES_ONLY */

#ifdef CK_TEMPLATES_ONLY
  #define _CHARMXI_CLASS_NAME EnzoReconstructorPLM<PLM_AthenaLimiter>
  PUPable_def_template(_CHARMXI_CLASS_NAME)
  #undef _CHARMXI_CLASS_NAME
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodStarMaker)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodStarMakerStochasticSF)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoMethodStarMakerSTARSS)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoBfieldMethodCT)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoObjectFeedbackSphere)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoPhysicsCosmology)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoPhysicsFluidProps)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoProblem)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoProlong)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoProlongMC1)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoProlongPoisson)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoRestrict)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoSolverCg)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoSolverDd)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoSolverDiagonal)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoSolverBiCgStab)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoSolverMg0)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoSolverJacobi)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoStopping)
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
  PUPable_def(EnzoUnits)
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: message EnzoMsgCheck;
 */
#ifndef CK_TEMPLATES_ONLY
void *CMessage_EnzoMsgCheck::operator new(size_t s){
  return EnzoMsgCheck::alloc(__idx, s, 0, 0, GroupDepNum{});
}
void *CMessage_EnzoMsgCheck::operator new(size_t s, int* sz){
  return EnzoMsgCheck::alloc(__idx, s, sz, 0, GroupDepNum{});
}
void *CMessage_EnzoMsgCheck::operator new(size_t s, int* sz,const int pb){
  return EnzoMsgCheck::alloc(__idx, s, sz, pb, GroupDepNum{});
}
void *CMessage_EnzoMsgCheck::operator new(size_t s, int* sz,const int pb, const GroupDepNum groupDepNum){
  return EnzoMsgCheck::alloc(__idx, s, sz, pb, groupDepNum);
}
void *CMessage_EnzoMsgCheck::operator new(size_t s, const int p) {
  return EnzoMsgCheck::alloc(__idx, s, 0, p, GroupDepNum{});
}
void *CMessage_EnzoMsgCheck::operator new(size_t s, const int p, const GroupDepNum groupDepNum) {
  return EnzoMsgCheck::alloc(__idx, s, 0, p, groupDepNum);
}
void* CMessage_EnzoMsgCheck::alloc(int msgnum, size_t sz, int *sizes, int pb, GroupDepNum groupDepNum) {
  CkpvAccess(_offsets)[0] = ALIGN_DEFAULT(sz);
  return CkAllocMsg(msgnum, CkpvAccess(_offsets)[0], pb, groupDepNum);
}
CMessage_EnzoMsgCheck::CMessage_EnzoMsgCheck() {
EnzoMsgCheck *newmsg = (EnzoMsgCheck *)this;
}
void CMessage_EnzoMsgCheck::dealloc(void *p) {
  CkFreeMsg(p);
}
void* CMessage_EnzoMsgCheck::pack(EnzoMsgCheck *msg) {
  return (void *) msg;
}
EnzoMsgCheck* CMessage_EnzoMsgCheck::unpack(void* buf) {
  EnzoMsgCheck *msg = (EnzoMsgCheck *) buf;
  return msg;
}
int CMessage_EnzoMsgCheck::__idx=0;
#endif /* CK_TEMPLATES_ONLY */


/* DEFS: readonly CProxy_EnzoSimulation proxy_enzo_simulation;
 */
extern CProxy_EnzoSimulation proxy_enzo_simulation;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_proxy_enzo_simulation(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|proxy_enzo_simulation;
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly CProxy_IoEnzoReader proxy_io_enzo_reader;
 */
extern CProxy_IoEnzoReader proxy_io_enzo_reader;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_proxy_io_enzo_reader(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|proxy_io_enzo_reader;
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly CProxy_IoEnzoWriter proxy_io_enzo_writer;
 */
extern CProxy_IoEnzoWriter proxy_io_enzo_writer;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_proxy_io_enzo_writer(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|proxy_io_enzo_writer;
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: readonly CProxy_EnzoLevelArray proxy_level_array;
 */
extern CProxy_EnzoLevelArray proxy_level_array;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_proxy_level_array(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|proxy_level_array;
}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: group EnzoSimulation: Simulation{
EnzoSimulation(const char *filename, int n);
void r_startup_begun(CkReductionMsg* impl_msg);
void p_get_msg_refine(const Index &index);
void p_get_msg_check(const Index &index);
void r_method_balance_count(CkReductionMsg* impl_msg);
void p_method_balance_check();
void r_method_check_enter(CkReductionMsg* impl_msg);
void p_check_done();
void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
void p_infer_set_array_count(int count);
void p_infer_array_created();
void p_infer_done();
void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
void p_io_reader_created();
void p_restart_next_level();
void p_restart_level_created();
void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
void p_fbnet_concatenate_sphere_lists();
void p_fbnet_done();
EnzoSimulation(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_EnzoSimulation::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoSimulation(const char *filename, int n);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_startup_begun(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoSimulation::r_startup_begun(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_get_msg_refine(const Index &index);
 */
void CProxyElement_EnzoSimulation::p_get_msg_refine(const Index &index, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_get_msg_check(const Index &index);
 */
void CProxyElement_EnzoSimulation::p_get_msg_check(const Index &index, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_balance_count(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoSimulation::r_method_balance_count(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_balance_check();
 */
void CProxyElement_EnzoSimulation::p_method_balance_check(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_check_enter(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoSimulation::r_method_check_enter(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_done();
 */
void CProxyElement_EnzoSimulation::p_check_done(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
 */
void CProxyElement_EnzoSimulation::p_set_io_writer(const CProxy_IoEnzoWriter &proxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_IoEnzoWriter &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_infer_set_array_count(int count);
 */
void CProxyElement_EnzoSimulation::p_infer_set_array_count(int count, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_infer_array_created();
 */
void CProxyElement_EnzoSimulation::p_infer_array_created(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_infer_done();
 */
void CProxyElement_EnzoSimulation::p_infer_done(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
 */
void CProxyElement_EnzoSimulation::p_set_io_reader(const CProxy_IoEnzoReader &proxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_IoEnzoReader &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_io_reader_created();
 */
void CProxyElement_EnzoSimulation::p_io_reader_created(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_next_level();
 */
void CProxyElement_EnzoSimulation::p_restart_next_level(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_level_created();
 */
void CProxyElement_EnzoSimulation::p_restart_level_created(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
 */
void CProxyElement_EnzoSimulation::p_set_level_array(const CProxy_EnzoLevelArray &proxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_EnzoLevelArray &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_fbnet_concatenate_sphere_lists();
 */
void CProxyElement_EnzoSimulation::p_fbnet_concatenate_sphere_lists(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_fbnet_done();
 */
void CProxyElement_EnzoSimulation::p_fbnet_done(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoSimulation(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoSimulation(const char *filename, int n);
 */
CkGroupID CProxy_EnzoSimulation::ckNew(const char *filename, int n, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const char *filename, int n
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_filename, impl_cnt_filename;
  impl_off_filename=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_filename=sizeof(char)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_filename;
    implP|impl_cnt_filename;
    implP|n;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_filename;
    implP|impl_cnt_filename;
    implP|n;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_filename,filename,impl_cnt_filename);
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  CkGroupID gId = CkCreateGroup(CkIndex_EnzoSimulation::__idx, CkIndex_EnzoSimulation::idx_EnzoSimulation_marshall1(), impl_msg);
  return gId;
}
  CProxy_EnzoSimulation::CProxy_EnzoSimulation(const char *filename, int n, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const char *filename, int n
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_filename, impl_cnt_filename;
  impl_off_filename=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_filename=sizeof(char)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_filename;
    implP|impl_cnt_filename;
    implP|n;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_filename;
    implP|impl_cnt_filename;
    implP|n;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_filename,filename,impl_cnt_filename);
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  ckSetGroupID(CkCreateGroup(CkIndex_EnzoSimulation::__idx, CkIndex_EnzoSimulation::idx_EnzoSimulation_marshall1(), impl_msg));
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_EnzoSimulation_marshall1() {
  int epidx = CkRegisterEp("EnzoSimulation(const char *filename, int n)",
      reinterpret_cast<CkCallFnPtr>(_call_EnzoSimulation_marshall1), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_EnzoSimulation_marshall1);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_EnzoSimulation_marshall1);

  return epidx;
}

void CkIndex_EnzoSimulation::_call_EnzoSimulation_marshall1(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const char *filename, int n*/
  PUP::fromMem implP(impl_buf);
  int impl_off_filename, impl_cnt_filename;
  implP|impl_off_filename;
  implP|impl_cnt_filename;
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *filename=(char *)(impl_buf+impl_off_filename);
  new (impl_obj_void) EnzoSimulation(filename, std::move(n.t));
}
int CkIndex_EnzoSimulation::_callmarshall_EnzoSimulation_marshall1(char* impl_buf, void* impl_obj_void) {
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const char *filename, int n*/
  PUP::fromMem implP(impl_buf);
  int impl_off_filename, impl_cnt_filename;
  implP|impl_off_filename;
  implP|impl_cnt_filename;
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *filename=(char *)(impl_buf+impl_off_filename);
  new (impl_obj_void) EnzoSimulation(filename, std::move(n.t));
  return implP.size();
}
void CkIndex_EnzoSimulation::_marshallmessagepup_EnzoSimulation_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const char *filename, int n*/
  PUP::fromMem implP(impl_buf);
  int impl_off_filename, impl_cnt_filename;
  implP|impl_off_filename;
  implP|impl_cnt_filename;
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *filename=(char *)(impl_buf+impl_off_filename);
  if (implDestP.hasComments()) implDestP.comment("filename");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*filename))<impl_cnt_filename;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|filename[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("n");
  implDestP|n;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_startup_begun(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoSimulation::r_startup_begun(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::r_startup_begun(CkReductionMsg* impl_msg, int npes, int *pes) {
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::r_startup_begun(CkReductionMsg* impl_msg, CmiGroup &grp) {
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_r_startup_begun_CkReductionMsg() {
  int epidx = CkRegisterEp("r_startup_begun(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_startup_begun_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_r_startup_begun_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->r_startup_begun((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_get_msg_refine(const Index &index);
 */
void CProxy_EnzoSimulation::p_get_msg_refine(const Index &index, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_get_msg_refine(const Index &index, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_get_msg_refine(const Index &index, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_get_msg_refine_marshall3() {
  int epidx = CkRegisterEp("p_get_msg_refine(const Index &index)",
      reinterpret_cast<CkCallFnPtr>(_call_p_get_msg_refine_marshall3), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_get_msg_refine_marshall3);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_get_msg_refine_marshall3);

  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_get_msg_refine_marshall3(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &index*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_get_msg_refine(std::move(index.t));
}
int CkIndex_EnzoSimulation::_callmarshall_p_get_msg_refine_marshall3(char* impl_buf, void* impl_obj_void) {
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const Index &index*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_get_msg_refine(std::move(index.t));
  return implP.size();
}
void CkIndex_EnzoSimulation::_marshallmessagepup_p_get_msg_refine_marshall3(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &index*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("index");
  implDestP|index;
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_get_msg_refine_3_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_get_msg_check(const Index &index);
 */
void CProxy_EnzoSimulation::p_get_msg_check(const Index &index, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_get_msg_check(const Index &index, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_get_msg_check(const Index &index, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_get_msg_check_marshall4() {
  int epidx = CkRegisterEp("p_get_msg_check(const Index &index)",
      reinterpret_cast<CkCallFnPtr>(_call_p_get_msg_check_marshall4), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_get_msg_check_marshall4);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_get_msg_check_marshall4);

  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_get_msg_check_marshall4(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &index*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_get_msg_check(std::move(index.t));
}
int CkIndex_EnzoSimulation::_callmarshall_p_get_msg_check_marshall4(char* impl_buf, void* impl_obj_void) {
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const Index &index*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_get_msg_check(std::move(index.t));
  return implP.size();
}
void CkIndex_EnzoSimulation::_marshallmessagepup_p_get_msg_check_marshall4(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &index*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> index;
  implP|index;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("index");
  implDestP|index;
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_get_msg_check_4_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_balance_count(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoSimulation::r_method_balance_count(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::r_method_balance_count(CkReductionMsg* impl_msg, int npes, int *pes) {
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::r_method_balance_count(CkReductionMsg* impl_msg, CmiGroup &grp) {
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_r_method_balance_count_CkReductionMsg() {
  int epidx = CkRegisterEp("r_method_balance_count(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_method_balance_count_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_r_method_balance_count_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->r_method_balance_count((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_balance_check();
 */
void CProxy_EnzoSimulation::p_method_balance_check(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_method_balance_check(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_method_balance_check(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_method_balance_check_void() {
  int epidx = CkRegisterEp("p_method_balance_check()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_balance_check_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_method_balance_check_void(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->p_method_balance_check();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_method_balance_check_6_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_check_enter(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoSimulation::r_method_check_enter(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::r_method_check_enter(CkReductionMsg* impl_msg, int npes, int *pes) {
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::r_method_check_enter(CkReductionMsg* impl_msg, CmiGroup &grp) {
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_r_method_check_enter_CkReductionMsg() {
  int epidx = CkRegisterEp("r_method_check_enter(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_method_check_enter_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_r_method_check_enter_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->r_method_check_enter((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_done();
 */
void CProxy_EnzoSimulation::p_check_done(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_check_done(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_check_done(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_check_done_void() {
  int epidx = CkRegisterEp("p_check_done()",
      reinterpret_cast<CkCallFnPtr>(_call_p_check_done_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_check_done_void(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->p_check_done();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_check_done_8_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
 */
void CProxy_EnzoSimulation::p_set_io_writer(const CProxy_IoEnzoWriter &proxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_IoEnzoWriter &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_set_io_writer(const CProxy_IoEnzoWriter &proxy, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: const CProxy_IoEnzoWriter &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
  }
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_set_io_writer(const CProxy_IoEnzoWriter &proxy, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: const CProxy_IoEnzoWriter &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
  }
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_set_io_writer_marshall9() {
  int epidx = CkRegisterEp("p_set_io_writer(const CProxy_IoEnzoWriter &proxy)",
      reinterpret_cast<CkCallFnPtr>(_call_p_set_io_writer_marshall9), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_set_io_writer_marshall9);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_set_io_writer_marshall9);

  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_set_io_writer_marshall9(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CProxy_IoEnzoWriter &proxy*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_IoEnzoWriter> proxy;
  implP|proxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_set_io_writer(std::move(proxy.t));
}
int CkIndex_EnzoSimulation::_callmarshall_p_set_io_writer_marshall9(char* impl_buf, void* impl_obj_void) {
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const CProxy_IoEnzoWriter &proxy*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_IoEnzoWriter> proxy;
  implP|proxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_set_io_writer(std::move(proxy.t));
  return implP.size();
}
void CkIndex_EnzoSimulation::_marshallmessagepup_p_set_io_writer_marshall9(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CProxy_IoEnzoWriter &proxy*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_IoEnzoWriter> proxy;
  implP|proxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("proxy");
  implDestP|proxy;
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_set_io_writer_9_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_infer_set_array_count(int count);
 */
void CProxy_EnzoSimulation::p_infer_set_array_count(int count, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_infer_set_array_count(int count, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count;
  }
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_infer_set_array_count(int count, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count;
  }
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_infer_set_array_count_marshall10() {
  int epidx = CkRegisterEp("p_infer_set_array_count(int count)",
      reinterpret_cast<CkCallFnPtr>(_call_p_infer_set_array_count_marshall10), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_infer_set_array_count_marshall10);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_infer_set_array_count_marshall10);

  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_infer_set_array_count_marshall10(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_infer_set_array_count(std::move(count.t));
}
int CkIndex_EnzoSimulation::_callmarshall_p_infer_set_array_count_marshall10(char* impl_buf, void* impl_obj_void) {
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_infer_set_array_count(std::move(count.t));
  return implP.size();
}
void CkIndex_EnzoSimulation::_marshallmessagepup_p_infer_set_array_count_marshall10(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("count");
  implDestP|count;
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_infer_set_array_count_10_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_infer_array_created();
 */
void CProxy_EnzoSimulation::p_infer_array_created(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_infer_array_created(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_infer_array_created(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_infer_array_created_void() {
  int epidx = CkRegisterEp("p_infer_array_created()",
      reinterpret_cast<CkCallFnPtr>(_call_p_infer_array_created_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_infer_array_created_void(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->p_infer_array_created();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_infer_array_created_11_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_infer_done();
 */
void CProxy_EnzoSimulation::p_infer_done(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_infer_done(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_infer_done(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_infer_done_void() {
  int epidx = CkRegisterEp("p_infer_done()",
      reinterpret_cast<CkCallFnPtr>(_call_p_infer_done_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_infer_done_void(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->p_infer_done();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_infer_done_12_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
 */
void CProxy_EnzoSimulation::p_set_io_reader(const CProxy_IoEnzoReader &proxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_IoEnzoReader &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_set_io_reader(const CProxy_IoEnzoReader &proxy, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: const CProxy_IoEnzoReader &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
  }
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_set_io_reader(const CProxy_IoEnzoReader &proxy, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: const CProxy_IoEnzoReader &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
  }
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_set_io_reader_marshall13() {
  int epidx = CkRegisterEp("p_set_io_reader(const CProxy_IoEnzoReader &proxy)",
      reinterpret_cast<CkCallFnPtr>(_call_p_set_io_reader_marshall13), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_set_io_reader_marshall13);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_set_io_reader_marshall13);

  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_set_io_reader_marshall13(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CProxy_IoEnzoReader &proxy*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_IoEnzoReader> proxy;
  implP|proxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_set_io_reader(std::move(proxy.t));
}
int CkIndex_EnzoSimulation::_callmarshall_p_set_io_reader_marshall13(char* impl_buf, void* impl_obj_void) {
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const CProxy_IoEnzoReader &proxy*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_IoEnzoReader> proxy;
  implP|proxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_set_io_reader(std::move(proxy.t));
  return implP.size();
}
void CkIndex_EnzoSimulation::_marshallmessagepup_p_set_io_reader_marshall13(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CProxy_IoEnzoReader &proxy*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_IoEnzoReader> proxy;
  implP|proxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("proxy");
  implDestP|proxy;
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_set_io_reader_13_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_io_reader_created();
 */
void CProxy_EnzoSimulation::p_io_reader_created(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_io_reader_created(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_io_reader_created(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_io_reader_created_void() {
  int epidx = CkRegisterEp("p_io_reader_created()",
      reinterpret_cast<CkCallFnPtr>(_call_p_io_reader_created_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_io_reader_created_void(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->p_io_reader_created();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_io_reader_created_14_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_next_level();
 */
void CProxy_EnzoSimulation::p_restart_next_level(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_restart_next_level(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_restart_next_level(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_restart_next_level_void() {
  int epidx = CkRegisterEp("p_restart_next_level()",
      reinterpret_cast<CkCallFnPtr>(_call_p_restart_next_level_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_restart_next_level_void(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->p_restart_next_level();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_restart_next_level_15_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_level_created();
 */
void CProxy_EnzoSimulation::p_restart_level_created(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_restart_level_created(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_restart_level_created(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_restart_level_created_void() {
  int epidx = CkRegisterEp("p_restart_level_created()",
      reinterpret_cast<CkCallFnPtr>(_call_p_restart_level_created_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_restart_level_created_void(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->p_restart_level_created();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_restart_level_created_16_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
 */
void CProxy_EnzoSimulation::p_set_level_array(const CProxy_EnzoLevelArray &proxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_EnzoLevelArray &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_set_level_array(const CProxy_EnzoLevelArray &proxy, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: const CProxy_EnzoLevelArray &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
  }
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_set_level_array(const CProxy_EnzoLevelArray &proxy, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: const CProxy_EnzoLevelArray &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
  }
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_set_level_array_marshall17() {
  int epidx = CkRegisterEp("p_set_level_array(const CProxy_EnzoLevelArray &proxy)",
      reinterpret_cast<CkCallFnPtr>(_call_p_set_level_array_marshall17), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_set_level_array_marshall17);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_set_level_array_marshall17);

  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_set_level_array_marshall17(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CProxy_EnzoLevelArray &proxy*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_EnzoLevelArray> proxy;
  implP|proxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_set_level_array(std::move(proxy.t));
}
int CkIndex_EnzoSimulation::_callmarshall_p_set_level_array_marshall17(char* impl_buf, void* impl_obj_void) {
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const CProxy_EnzoLevelArray &proxy*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_EnzoLevelArray> proxy;
  implP|proxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_set_level_array(std::move(proxy.t));
  return implP.size();
}
void CkIndex_EnzoSimulation::_marshallmessagepup_p_set_level_array_marshall17(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CProxy_EnzoLevelArray &proxy*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_EnzoLevelArray> proxy;
  implP|proxy;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("proxy");
  implDestP|proxy;
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_set_level_array_17_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_fbnet_concatenate_sphere_lists();
 */
void CProxy_EnzoSimulation::p_fbnet_concatenate_sphere_lists(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_fbnet_concatenate_sphere_lists(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_fbnet_concatenate_sphere_lists(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_fbnet_concatenate_sphere_lists_void() {
  int epidx = CkRegisterEp("p_fbnet_concatenate_sphere_lists()",
      reinterpret_cast<CkCallFnPtr>(_call_p_fbnet_concatenate_sphere_lists_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_fbnet_concatenate_sphere_lists_void(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->p_fbnet_concatenate_sphere_lists();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_fbnet_concatenate_sphere_lists_18_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_fbnet_done();
 */
void CProxy_EnzoSimulation::p_fbnet_done(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_EnzoSimulation::p_fbnet_done(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_EnzoSimulation::p_fbnet_done(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_EnzoSimulation::reg_p_fbnet_done_void() {
  int epidx = CkRegisterEp("p_fbnet_done()",
      reinterpret_cast<CkCallFnPtr>(_call_p_fbnet_done_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_p_fbnet_done_void(void* impl_msg, void* impl_obj_void)
{
  EnzoSimulation* impl_obj = static_cast<EnzoSimulation*>(impl_obj_void);
  impl_obj->p_fbnet_done();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoSimulation::p_fbnet_done_19_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoSimulation(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_EnzoSimulation::reg_EnzoSimulation_CkMigrateMessage() {
  int epidx = CkRegisterEp("EnzoSimulation(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_EnzoSimulation_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoSimulation::_call_EnzoSimulation_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  new (impl_obj_void) EnzoSimulation((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoSimulation(const char *filename, int n);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_startup_begun(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoSimulation::r_startup_begun(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_r_startup_begun_CkReductionMsg(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_get_msg_refine(const Index &index);
 */
void CProxySection_EnzoSimulation::p_get_msg_refine(const Index &index, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_get_msg_refine_marshall3(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_get_msg_check(const Index &index);
 */
void CProxySection_EnzoSimulation::p_get_msg_check(const Index &index, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const Index &index
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)index;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_get_msg_check_marshall4(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_balance_count(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoSimulation::r_method_balance_count(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_r_method_balance_count_CkReductionMsg(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_balance_check();
 */
void CProxySection_EnzoSimulation::p_method_balance_check(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_method_balance_check_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_check_enter(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoSimulation::r_method_check_enter(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_r_method_check_enter_CkReductionMsg(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_done();
 */
void CProxySection_EnzoSimulation::p_check_done(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_check_done_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
 */
void CProxySection_EnzoSimulation::p_set_io_writer(const CProxy_IoEnzoWriter &proxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_IoEnzoWriter &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoWriter>::type>::type &)proxy;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_set_io_writer_marshall9(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_infer_set_array_count(int count);
 */
void CProxySection_EnzoSimulation::p_infer_set_array_count(int count, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_infer_set_array_count_marshall10(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_infer_array_created();
 */
void CProxySection_EnzoSimulation::p_infer_array_created(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_infer_array_created_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_infer_done();
 */
void CProxySection_EnzoSimulation::p_infer_done(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_infer_done_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
 */
void CProxySection_EnzoSimulation::p_set_io_reader(const CProxy_IoEnzoReader &proxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_IoEnzoReader &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_IoEnzoReader>::type>::type &)proxy;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_set_io_reader_marshall13(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_io_reader_created();
 */
void CProxySection_EnzoSimulation::p_io_reader_created(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_io_reader_created_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_next_level();
 */
void CProxySection_EnzoSimulation::p_restart_next_level(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_restart_next_level_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_level_created();
 */
void CProxySection_EnzoSimulation::p_restart_level_created(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_restart_level_created_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
 */
void CProxySection_EnzoSimulation::p_set_level_array(const CProxy_EnzoLevelArray &proxy, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_EnzoLevelArray &proxy
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_EnzoLevelArray>::type>::type &)proxy;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_set_level_array_marshall17(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_fbnet_concatenate_sphere_lists();
 */
void CProxySection_EnzoSimulation::p_fbnet_concatenate_sphere_lists(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_fbnet_concatenate_sphere_lists_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_fbnet_done();
 */
void CProxySection_EnzoSimulation::p_fbnet_done(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_EnzoSimulation::idx_p_fbnet_done_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoSimulation(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_EnzoSimulation::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_Simulation::__idx);
   CkRegisterGroupIrr(__idx,EnzoSimulation::isIrreducible());
  // REG: EnzoSimulation(const char *filename, int n);
  idx_EnzoSimulation_marshall1();

  // REG: void r_startup_begun(CkReductionMsg* impl_msg);
  idx_r_startup_begun_CkReductionMsg();

  // REG: void p_get_msg_refine(const Index &index);
  idx_p_get_msg_refine_marshall3();

  // REG: void p_get_msg_check(const Index &index);
  idx_p_get_msg_check_marshall4();

  // REG: void r_method_balance_count(CkReductionMsg* impl_msg);
  idx_r_method_balance_count_CkReductionMsg();

  // REG: void p_method_balance_check();
  idx_p_method_balance_check_void();

  // REG: void r_method_check_enter(CkReductionMsg* impl_msg);
  idx_r_method_check_enter_CkReductionMsg();

  // REG: void p_check_done();
  idx_p_check_done_void();

  // REG: void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
  idx_p_set_io_writer_marshall9();

  // REG: void p_infer_set_array_count(int count);
  idx_p_infer_set_array_count_marshall10();

  // REG: void p_infer_array_created();
  idx_p_infer_array_created_void();

  // REG: void p_infer_done();
  idx_p_infer_done_void();

  // REG: void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
  idx_p_set_io_reader_marshall13();

  // REG: void p_io_reader_created();
  idx_p_io_reader_created_void();

  // REG: void p_restart_next_level();
  idx_p_restart_next_level_void();

  // REG: void p_restart_level_created();
  idx_p_restart_level_created_void();

  // REG: void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
  idx_p_set_level_array_marshall17();

  // REG: void p_fbnet_concatenate_sphere_lists();
  idx_p_fbnet_concatenate_sphere_lists_void();

  // REG: void p_fbnet_done();
  idx_p_fbnet_done_void();

  // REG: EnzoSimulation(CkMigrateMessage* impl_msg);
  idx_EnzoSimulation_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_EnzoSimulation_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: array EnzoBlock: Block{
EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
void p_set_msg_refine(MsgRefine* impl_msg);
void p_set_msg_check(EnzoMsgCheck* impl_msg);
EnzoBlock();
void p_method_feedback_starss_end();
void p_method_m1_closure_solve_transport_eqn();
void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
void r_method_turbulence_end(CkReductionMsg* impl_msg);
void p_initial_hdf5_recv(MsgInitial* impl_msg);
void p_method_balance_migrate();
void p_method_balance_done();
void p_method_gravity_continue();
void p_method_gravity_end();
void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
void p_method_infer_count_arrays(int count);
void p_method_infer_request_data(const int *il3);
void p_method_infer_update(int n, const char *buffer, const int *il3);
void p_method_infer_exit();
void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
void p_method_fbnet_exit();
void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
void p_check_write_next(int num_files, const std::string &ordering);
void p_check_done();
void p_restart_refine(const int *ic3, int io_reader);
void p_restart_set_data(EnzoMsgCheck* impl_msg);
void p_restart_done();
void p_method_accretion_end();
void p_solver_cg_matvec();
void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
void p_solver_cg_loop_2();
void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
void p_solver_bicgstab_loop_2();
void p_solver_bicgstab_loop_3();
void p_solver_bicgstab_loop_8();
void p_solver_bicgstab_loop_9();
void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
void p_solver_dd_solve_coarse();
void p_solver_dd_solve_domain();
void p_solver_dd_last_smooth();
void r_solver_dd_barrier(CkReductionMsg* impl_msg);
void r_solver_dd_end(CkReductionMsg* impl_msg);
void p_solver_jacobi_continue();
void p_solver_mg0_restrict();
void p_solver_mg0_solve_coarse();
void p_solver_mg0_post_smooth();
void p_solver_mg0_last_smooth();
void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
EnzoBlock(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_EnzoBlock::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CProxySection_EnzoBlock::contribute(CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, userData, fragSize);
}

void CProxySection_EnzoBlock::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, userData, fragSize);
}

template <typename T>
void CProxySection_EnzoBlock::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, userData, fragSize);
}

void CProxySection_EnzoBlock::contribute(CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, cb, userData, fragSize);
}

void CProxySection_EnzoBlock::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, cb, userData, fragSize);
}

template <typename T>
void CProxySection_EnzoBlock::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, cb, userData, fragSize);
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
 */
void CProxyElement_EnzoBlock::insert(const process_type &ip_source, const MsgType &msg_type, int onPE, const CkEntryOptions *impl_e_opts)
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
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_EnzoBlock::idx_EnzoBlock_marshall1(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
void CProxyElement_EnzoBlock::p_set_msg_refine(MsgRefine* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_set_msg_refine_MsgRefine(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_msg_check(EnzoMsgCheck* impl_msg);
 */
void CProxyElement_EnzoBlock::p_set_msg_check(EnzoMsgCheck* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_set_msg_check_EnzoMsgCheck(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoBlock();
 */
void CProxyElement_EnzoBlock::insert(int onPE, const CkEntryOptions *impl_e_opts)
{ 
   void *impl_msg = CkAllocSysMsg(impl_e_opts);
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_EnzoBlock::idx_EnzoBlock_void(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_feedback_starss_end();
 */
void CProxyElement_EnzoBlock::p_method_feedback_starss_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_feedback_starss_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_m1_closure_solve_transport_eqn();
 */
void CProxyElement_EnzoBlock::p_method_m1_closure_solve_transport_eqn(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_m1_closure_solve_transport_eqn_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_m1_closure_set_global_averages_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_turbulence_end(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_method_turbulence_end(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_method_turbulence_end_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_initial_hdf5_recv(MsgInitial* impl_msg);
 */
void CProxyElement_EnzoBlock::p_initial_hdf5_recv(MsgInitial* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_initial_hdf5_recv_MsgInitial(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_balance_migrate();
 */
void CProxyElement_EnzoBlock::p_method_balance_migrate(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_balance_migrate_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_balance_done();
 */
void CProxyElement_EnzoBlock::p_method_balance_done(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_balance_done_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_gravity_continue();
 */
void CProxyElement_EnzoBlock::p_method_gravity_continue(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_gravity_continue_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_gravity_end();
 */
void CProxyElement_EnzoBlock::p_method_gravity_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_gravity_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
 */
void CProxyElement_EnzoBlock::p_method_infer_merge_masks(int n, const char *mask, const int *ic3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const char *mask, const int *ic3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_mask, impl_cnt_mask;
  impl_off_mask=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_mask=sizeof(char)*(n));
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_mask;
    implP|impl_cnt_mask;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_mask;
    implP|impl_cnt_mask;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_mask,mask,impl_cnt_mask);
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_merge_masks_marshall14(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_count_arrays(int count);
 */
void CProxyElement_EnzoBlock::p_method_infer_count_arrays(int count, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_count_arrays_marshall15(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_request_data(const int *il3);
 */
void CProxyElement_EnzoBlock::p_method_infer_request_data(const int *il3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const int *il3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_il3, impl_cnt_il3;
  impl_off_il3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_il3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_il3;
    implP|impl_cnt_il3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_il3;
    implP|impl_cnt_il3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_il3,il3,impl_cnt_il3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_request_data_marshall16(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_update(int n, const char *buffer, const int *il3);
 */
void CProxyElement_EnzoBlock::p_method_infer_update(int n, const char *buffer, const int *il3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const char *buffer, const int *il3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_buffer, impl_cnt_buffer;
  impl_off_buffer=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_buffer=sizeof(char)*(n));
  int impl_off_il3, impl_cnt_il3;
  impl_off_il3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_il3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    implP|impl_off_il3;
    implP|impl_cnt_il3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    implP|impl_off_il3;
    implP|impl_cnt_il3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_buffer,buffer,impl_cnt_buffer);
  memcpy(impl_buf+impl_off_il3,il3,impl_cnt_il3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_update_marshall17(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_exit();
 */
void CProxyElement_EnzoBlock::p_method_infer_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::p_method_fbnet_update_mesh(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_fbnet_update_mesh_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_fbnet_exit();
 */
void CProxyElement_EnzoBlock::p_method_fbnet_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_fbnet_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
 */
void CProxyElement_EnzoBlock::p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int num_files, const std::string &ordering, const std::string &name_dir
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)name_dir;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)name_dir;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_check_write_first_marshall21(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_write_next(int num_files, const std::string &ordering);
 */
void CProxyElement_EnzoBlock::p_check_write_next(int num_files, const std::string &ordering, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int num_files, const std::string &ordering
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_check_write_next_marshall22(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_done();
 */
void CProxyElement_EnzoBlock::p_check_done(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_check_done_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_refine(const int *ic3, int io_reader);
 */
void CProxyElement_EnzoBlock::p_restart_refine(const int *ic3, int io_reader, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const int *ic3, int io_reader
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|io_reader;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|io_reader;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_restart_refine_marshall24(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_set_data(EnzoMsgCheck* impl_msg);
 */
void CProxyElement_EnzoBlock::p_restart_set_data(EnzoMsgCheck* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_restart_set_data_EnzoMsgCheck(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_done();
 */
void CProxyElement_EnzoBlock::p_restart_done(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_restart_done_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_accretion_end();
 */
void CProxyElement_EnzoBlock::p_method_accretion_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_accretion_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_cg_matvec();
 */
void CProxyElement_EnzoBlock::p_solver_cg_matvec(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_cg_matvec_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_cg_loop_0a(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_0a_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_cg_loop_0b(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_0b_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_cg_shift_1(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_shift_1_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_cg_loop_2();
 */
void CProxyElement_EnzoBlock::p_solver_cg_loop_2(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_cg_loop_2_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_cg_loop_3(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_3_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_cg_loop_5(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_5_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_bicgstab_start_1(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_start_1_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_bicgstab_start_3(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_start_3_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_5_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_11_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_13_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_15_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_2();
 */
void CProxyElement_EnzoBlock::p_solver_bicgstab_loop_2(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_2_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_3();
 */
void CProxyElement_EnzoBlock::p_solver_bicgstab_loop_3(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_3_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_8();
 */
void CProxyElement_EnzoBlock::p_solver_bicgstab_loop_8(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_8_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_9();
 */
void CProxyElement_EnzoBlock::p_solver_bicgstab_loop_9(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_9_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
 */
void CProxyElement_EnzoBlock::p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_dot, impl_cnt_dot;
  impl_off_dot=impl_off=CK_ALIGN(impl_off,sizeof(long double));
  impl_off+=(impl_cnt_dot=sizeof(long double)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
    implP|iter;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
    implP|iter;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_dot,dot,impl_cnt_dot);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_dot_recv_parent_marshall45(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
 */
void CProxyElement_EnzoBlock::p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const long double *dot, const std::vector<int> &isa, int i_function
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_dot, impl_cnt_dot;
  impl_off_dot=impl_off=CK_ALIGN(impl_off,sizeof(long double));
  impl_off+=(impl_cnt_dot=sizeof(long double)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_dot,dot,impl_cnt_dot);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_dot_recv_children_marshall46(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::p_solver_dd_restrict_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_restrict_recv_FieldMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::p_solver_dd_prolong_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_prolong_recv_FieldMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_solve_coarse();
 */
void CProxyElement_EnzoBlock::p_solver_dd_solve_coarse(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_solve_coarse_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_solve_domain();
 */
void CProxyElement_EnzoBlock::p_solver_dd_solve_domain(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_solve_domain_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_last_smooth();
 */
void CProxyElement_EnzoBlock::p_solver_dd_last_smooth(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_last_smooth_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_dd_barrier(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_dd_barrier(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_dd_barrier_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_dd_end(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_dd_end(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_dd_end_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_jacobi_continue();
 */
void CProxyElement_EnzoBlock::p_solver_jacobi_continue(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_jacobi_continue_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_restrict();
 */
void CProxyElement_EnzoBlock::p_solver_mg0_restrict(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_restrict_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_solve_coarse();
 */
void CProxyElement_EnzoBlock::p_solver_mg0_solve_coarse(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_solve_coarse_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_post_smooth();
 */
void CProxyElement_EnzoBlock::p_solver_mg0_post_smooth(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_post_smooth_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_last_smooth();
 */
void CProxyElement_EnzoBlock::p_solver_mg0_last_smooth(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_last_smooth_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_mg0_begin_solve(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_mg0_begin_solve_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::r_solver_mg0_barrier(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_mg0_barrier_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::p_solver_mg0_prolong_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_prolong_recv_FieldMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
 */
void CProxyElement_EnzoBlock::p_solver_mg0_restrict_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_restrict_recv_FieldMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoBlock(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
 */
CkArrayID CProxy_EnzoBlock::ckNew(const process_type &ip_source, const MsgType &msg_type, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
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
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_EnzoBlock::idx_EnzoBlock_marshall1(), opts);
  return gId;
}
void CProxy_EnzoBlock::ckNew(const process_type &ip_source, const MsgType &msg_type, const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
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
  CkSendAsyncCreateArray(CkIndex_EnzoBlock::idx_EnzoBlock_marshall1(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_EnzoBlock_marshall1() {
  int epidx = CkRegisterEp("EnzoBlock(const process_type &ip_source, const MsgType &msg_type)",
      reinterpret_cast<CkCallFnPtr>(_call_EnzoBlock_marshall1), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_EnzoBlock_marshall1);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_EnzoBlock_marshall1);

  return epidx;
}

void CkIndex_EnzoBlock::_call_EnzoBlock_marshall1(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
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
  new (impl_obj_void) EnzoBlock(std::move(ip_source.t), std::move(msg_type.t));
}
int CkIndex_EnzoBlock::_callmarshall_EnzoBlock_marshall1(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const process_type &ip_source, const MsgType &msg_type*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<process_type> ip_source;
  implP|ip_source;
  PUP::detail::TemporaryObjectHolder<MsgType> msg_type;
  implP|msg_type;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) EnzoBlock(std::move(ip_source.t), std::move(msg_type.t));
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_EnzoBlock_marshall1(PUP::er &implDestP,void *impl_msg) {
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
void CProxy_EnzoBlock::p_set_msg_refine(MsgRefine* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_set_msg_refine_MsgRefine(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_set_msg_refine_MsgRefine() {
  int epidx = CkRegisterEp("p_set_msg_refine(MsgRefine* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_set_msg_refine_MsgRefine), CMessage_MsgRefine::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)MsgRefine::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_set_msg_refine_MsgRefine(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_set_msg_refine((MsgRefine*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_msg_check(EnzoMsgCheck* impl_msg);
 */
void CProxy_EnzoBlock::p_set_msg_check(EnzoMsgCheck* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_set_msg_check_EnzoMsgCheck(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_set_msg_check_EnzoMsgCheck() {
  int epidx = CkRegisterEp("p_set_msg_check(EnzoMsgCheck* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_set_msg_check_EnzoMsgCheck), CMessage_EnzoMsgCheck::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)EnzoMsgCheck::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_set_msg_check_EnzoMsgCheck(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_set_msg_check((EnzoMsgCheck*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoBlock();
 */
CkArrayID CProxy_EnzoBlock::ckNew(const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_EnzoBlock::idx_EnzoBlock_void(), opts);
  return gId;
}
void CProxy_EnzoBlock::ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_EnzoBlock::idx_EnzoBlock_void(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_EnzoBlock_void() {
  int epidx = CkRegisterEp("EnzoBlock()",
      reinterpret_cast<CkCallFnPtr>(_call_EnzoBlock_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_EnzoBlock_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  new (impl_obj_void) EnzoBlock();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_feedback_starss_end();
 */
void CProxy_EnzoBlock::p_method_feedback_starss_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_feedback_starss_end_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_feedback_starss_end_void() {
  int epidx = CkRegisterEp("p_method_feedback_starss_end()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_feedback_starss_end_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_feedback_starss_end_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_feedback_starss_end();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_feedback_starss_end_5_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_m1_closure_solve_transport_eqn();
 */
void CProxy_EnzoBlock::p_method_m1_closure_solve_transport_eqn(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_m1_closure_solve_transport_eqn_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_m1_closure_solve_transport_eqn_void() {
  int epidx = CkRegisterEp("p_method_m1_closure_solve_transport_eqn()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_m1_closure_solve_transport_eqn_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_m1_closure_solve_transport_eqn_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_m1_closure_solve_transport_eqn();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_m1_closure_solve_transport_eqn_6_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_m1_closure_set_global_averages_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_m1_closure_set_global_averages_CkReductionMsg() {
  int epidx = CkRegisterEp("p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_m1_closure_set_global_averages_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_m1_closure_set_global_averages_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_m1_closure_set_global_averages((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_turbulence_end(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_method_turbulence_end(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_method_turbulence_end_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_method_turbulence_end_CkReductionMsg() {
  int epidx = CkRegisterEp("r_method_turbulence_end(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_method_turbulence_end_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_method_turbulence_end_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_method_turbulence_end((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_initial_hdf5_recv(MsgInitial* impl_msg);
 */
void CProxy_EnzoBlock::p_initial_hdf5_recv(MsgInitial* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_initial_hdf5_recv_MsgInitial(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_initial_hdf5_recv_MsgInitial() {
  int epidx = CkRegisterEp("p_initial_hdf5_recv(MsgInitial* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_initial_hdf5_recv_MsgInitial), CMessage_MsgInitial::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)MsgInitial::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_initial_hdf5_recv_MsgInitial(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_initial_hdf5_recv((MsgInitial*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_balance_migrate();
 */
void CProxy_EnzoBlock::p_method_balance_migrate(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_balance_migrate_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_balance_migrate_void() {
  int epidx = CkRegisterEp("p_method_balance_migrate()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_balance_migrate_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_balance_migrate_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_balance_migrate();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_balance_migrate_10_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_balance_done();
 */
void CProxy_EnzoBlock::p_method_balance_done(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_balance_done_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_balance_done_void() {
  int epidx = CkRegisterEp("p_method_balance_done()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_balance_done_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_balance_done_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_balance_done();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_balance_done_11_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_gravity_continue();
 */
void CProxy_EnzoBlock::p_method_gravity_continue(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_gravity_continue_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_gravity_continue_void() {
  int epidx = CkRegisterEp("p_method_gravity_continue()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_gravity_continue_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_gravity_continue_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_gravity_continue();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_gravity_continue_12_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_gravity_end();
 */
void CProxy_EnzoBlock::p_method_gravity_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_gravity_end_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_gravity_end_void() {
  int epidx = CkRegisterEp("p_method_gravity_end()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_gravity_end_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_gravity_end_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_gravity_end();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_gravity_end_13_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
 */
void CProxy_EnzoBlock::p_method_infer_merge_masks(int n, const char *mask, const int *ic3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const char *mask, const int *ic3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_mask, impl_cnt_mask;
  impl_off_mask=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_mask=sizeof(char)*(n));
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_mask;
    implP|impl_cnt_mask;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_mask;
    implP|impl_cnt_mask;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_mask,mask,impl_cnt_mask);
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_merge_masks_marshall14(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_infer_merge_masks_marshall14() {
  int epidx = CkRegisterEp("p_method_infer_merge_masks(int n, const char *mask, const int *ic3)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_infer_merge_masks_marshall14), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_method_infer_merge_masks_marshall14);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_method_infer_merge_masks_marshall14);

  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_infer_merge_masks_marshall14(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const char *mask, const int *ic3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_mask, impl_cnt_mask;
  implP|impl_off_mask;
  implP|impl_cnt_mask;
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *mask=(char *)(impl_buf+impl_off_mask);
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  impl_obj->p_method_infer_merge_masks(std::move(n.t), mask, ic3);
}
int CkIndex_EnzoBlock::_callmarshall_p_method_infer_merge_masks_marshall14(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int n, const char *mask, const int *ic3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_mask, impl_cnt_mask;
  implP|impl_off_mask;
  implP|impl_cnt_mask;
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *mask=(char *)(impl_buf+impl_off_mask);
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  impl_obj->p_method_infer_merge_masks(std::move(n.t), mask, ic3);
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_p_method_infer_merge_masks_marshall14(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const char *mask, const int *ic3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_mask, impl_cnt_mask;
  implP|impl_off_mask;
  implP|impl_cnt_mask;
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *mask=(char *)(impl_buf+impl_off_mask);
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  if (implDestP.hasComments()) implDestP.comment("n");
  implDestP|n;
  if (implDestP.hasComments()) implDestP.comment("mask");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*mask))<impl_cnt_mask;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|mask[impl_i];
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
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_infer_merge_masks_14_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_count_arrays(int count);
 */
void CProxy_EnzoBlock::p_method_infer_count_arrays(int count, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_count_arrays_marshall15(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_infer_count_arrays_marshall15() {
  int epidx = CkRegisterEp("p_method_infer_count_arrays(int count)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_infer_count_arrays_marshall15), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_method_infer_count_arrays_marshall15);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_method_infer_count_arrays_marshall15);

  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_infer_count_arrays_marshall15(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_method_infer_count_arrays(std::move(count.t));
}
int CkIndex_EnzoBlock::_callmarshall_p_method_infer_count_arrays_marshall15(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_method_infer_count_arrays(std::move(count.t));
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_p_method_infer_count_arrays_marshall15(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int count*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> count;
  implP|count;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("count");
  implDestP|count;
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_infer_count_arrays_15_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_request_data(const int *il3);
 */
void CProxy_EnzoBlock::p_method_infer_request_data(const int *il3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const int *il3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_il3, impl_cnt_il3;
  impl_off_il3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_il3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_il3;
    implP|impl_cnt_il3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_il3;
    implP|impl_cnt_il3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_il3,il3,impl_cnt_il3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_request_data_marshall16(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_infer_request_data_marshall16() {
  int epidx = CkRegisterEp("p_method_infer_request_data(const int *il3)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_infer_request_data_marshall16), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_method_infer_request_data_marshall16);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_method_infer_request_data_marshall16);

  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_infer_request_data_marshall16(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const int *il3*/
  PUP::fromMem implP(impl_buf);
  int impl_off_il3, impl_cnt_il3;
  implP|impl_off_il3;
  implP|impl_cnt_il3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  int *il3=(int *)(impl_buf+impl_off_il3);
  impl_obj->p_method_infer_request_data(il3);
}
int CkIndex_EnzoBlock::_callmarshall_p_method_infer_request_data_marshall16(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const int *il3*/
  PUP::fromMem implP(impl_buf);
  int impl_off_il3, impl_cnt_il3;
  implP|impl_off_il3;
  implP|impl_cnt_il3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  int *il3=(int *)(impl_buf+impl_off_il3);
  impl_obj->p_method_infer_request_data(il3);
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_p_method_infer_request_data_marshall16(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const int *il3*/
  PUP::fromMem implP(impl_buf);
  int impl_off_il3, impl_cnt_il3;
  implP|impl_off_il3;
  implP|impl_cnt_il3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  int *il3=(int *)(impl_buf+impl_off_il3);
  if (implDestP.hasComments()) implDestP.comment("il3");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*il3))<impl_cnt_il3;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|il3[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_infer_request_data_16_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_update(int n, const char *buffer, const int *il3);
 */
void CProxy_EnzoBlock::p_method_infer_update(int n, const char *buffer, const int *il3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const char *buffer, const int *il3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_buffer, impl_cnt_buffer;
  impl_off_buffer=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_buffer=sizeof(char)*(n));
  int impl_off_il3, impl_cnt_il3;
  impl_off_il3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_il3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    implP|impl_off_il3;
    implP|impl_cnt_il3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    implP|impl_off_il3;
    implP|impl_cnt_il3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_buffer,buffer,impl_cnt_buffer);
  memcpy(impl_buf+impl_off_il3,il3,impl_cnt_il3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_update_marshall17(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_infer_update_marshall17() {
  int epidx = CkRegisterEp("p_method_infer_update(int n, const char *buffer, const int *il3)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_infer_update_marshall17), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_method_infer_update_marshall17);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_method_infer_update_marshall17);

  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_infer_update_marshall17(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const char *buffer, const int *il3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_buffer, impl_cnt_buffer;
  implP|impl_off_buffer;
  implP|impl_cnt_buffer;
  int impl_off_il3, impl_cnt_il3;
  implP|impl_off_il3;
  implP|impl_cnt_il3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *buffer=(char *)(impl_buf+impl_off_buffer);
  int *il3=(int *)(impl_buf+impl_off_il3);
  impl_obj->p_method_infer_update(std::move(n.t), buffer, il3);
}
int CkIndex_EnzoBlock::_callmarshall_p_method_infer_update_marshall17(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int n, const char *buffer, const int *il3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_buffer, impl_cnt_buffer;
  implP|impl_off_buffer;
  implP|impl_cnt_buffer;
  int impl_off_il3, impl_cnt_il3;
  implP|impl_off_il3;
  implP|impl_cnt_il3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *buffer=(char *)(impl_buf+impl_off_buffer);
  int *il3=(int *)(impl_buf+impl_off_il3);
  impl_obj->p_method_infer_update(std::move(n.t), buffer, il3);
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_p_method_infer_update_marshall17(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const char *buffer, const int *il3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_buffer, impl_cnt_buffer;
  implP|impl_off_buffer;
  implP|impl_cnt_buffer;
  int impl_off_il3, impl_cnt_il3;
  implP|impl_off_il3;
  implP|impl_cnt_il3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *buffer=(char *)(impl_buf+impl_off_buffer);
  int *il3=(int *)(impl_buf+impl_off_il3);
  if (implDestP.hasComments()) implDestP.comment("n");
  implDestP|n;
  if (implDestP.hasComments()) implDestP.comment("buffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*buffer))<impl_cnt_buffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|buffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("il3");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*il3))<impl_cnt_il3;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|il3[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_infer_update_17_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_exit();
 */
void CProxy_EnzoBlock::p_method_infer_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_exit_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_infer_exit_void() {
  int epidx = CkRegisterEp("p_method_infer_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_infer_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_infer_exit_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_infer_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_infer_exit_18_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::p_method_fbnet_update_mesh(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_fbnet_update_mesh_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_fbnet_update_mesh_CkReductionMsg() {
  int epidx = CkRegisterEp("p_method_fbnet_update_mesh(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_fbnet_update_mesh_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_fbnet_update_mesh_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_fbnet_update_mesh((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_fbnet_exit();
 */
void CProxy_EnzoBlock::p_method_fbnet_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_fbnet_exit_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_fbnet_exit_void() {
  int epidx = CkRegisterEp("p_method_fbnet_exit()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_fbnet_exit_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_fbnet_exit_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_fbnet_exit();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_fbnet_exit_20_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
 */
void CProxy_EnzoBlock::p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int num_files, const std::string &ordering, const std::string &name_dir
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)name_dir;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)name_dir;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_check_write_first_marshall21(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_check_write_first_marshall21() {
  int epidx = CkRegisterEp("p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir)",
      reinterpret_cast<CkCallFnPtr>(_call_p_check_write_first_marshall21), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_check_write_first_marshall21);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_check_write_first_marshall21);

  return epidx;
}

void CkIndex_EnzoBlock::_call_p_check_write_first_marshall21(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int num_files, const std::string &ordering, const std::string &name_dir*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> num_files;
  implP|num_files;
  PUP::detail::TemporaryObjectHolder<std::string> ordering;
  implP|ordering;
  PUP::detail::TemporaryObjectHolder<std::string> name_dir;
  implP|name_dir;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_check_write_first(std::move(num_files.t), std::move(ordering.t), std::move(name_dir.t));
}
int CkIndex_EnzoBlock::_callmarshall_p_check_write_first_marshall21(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int num_files, const std::string &ordering, const std::string &name_dir*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> num_files;
  implP|num_files;
  PUP::detail::TemporaryObjectHolder<std::string> ordering;
  implP|ordering;
  PUP::detail::TemporaryObjectHolder<std::string> name_dir;
  implP|name_dir;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_check_write_first(std::move(num_files.t), std::move(ordering.t), std::move(name_dir.t));
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_p_check_write_first_marshall21(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int num_files, const std::string &ordering, const std::string &name_dir*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> num_files;
  implP|num_files;
  PUP::detail::TemporaryObjectHolder<std::string> ordering;
  implP|ordering;
  PUP::detail::TemporaryObjectHolder<std::string> name_dir;
  implP|name_dir;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("num_files");
  implDestP|num_files;
  if (implDestP.hasComments()) implDestP.comment("ordering");
  implDestP|ordering;
  if (implDestP.hasComments()) implDestP.comment("name_dir");
  implDestP|name_dir;
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_check_write_first_21_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_write_next(int num_files, const std::string &ordering);
 */
void CProxy_EnzoBlock::p_check_write_next(int num_files, const std::string &ordering, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int num_files, const std::string &ordering
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_check_write_next_marshall22(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_check_write_next_marshall22() {
  int epidx = CkRegisterEp("p_check_write_next(int num_files, const std::string &ordering)",
      reinterpret_cast<CkCallFnPtr>(_call_p_check_write_next_marshall22), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_check_write_next_marshall22);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_check_write_next_marshall22);

  return epidx;
}

void CkIndex_EnzoBlock::_call_p_check_write_next_marshall22(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int num_files, const std::string &ordering*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> num_files;
  implP|num_files;
  PUP::detail::TemporaryObjectHolder<std::string> ordering;
  implP|ordering;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_check_write_next(std::move(num_files.t), std::move(ordering.t));
}
int CkIndex_EnzoBlock::_callmarshall_p_check_write_next_marshall22(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int num_files, const std::string &ordering*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> num_files;
  implP|num_files;
  PUP::detail::TemporaryObjectHolder<std::string> ordering;
  implP|ordering;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_check_write_next(std::move(num_files.t), std::move(ordering.t));
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_p_check_write_next_marshall22(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int num_files, const std::string &ordering*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> num_files;
  implP|num_files;
  PUP::detail::TemporaryObjectHolder<std::string> ordering;
  implP|ordering;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("num_files");
  implDestP|num_files;
  if (implDestP.hasComments()) implDestP.comment("ordering");
  implDestP|ordering;
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_check_write_next_22_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_done();
 */
void CProxy_EnzoBlock::p_check_done(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_check_done_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_check_done_void() {
  int epidx = CkRegisterEp("p_check_done()",
      reinterpret_cast<CkCallFnPtr>(_call_p_check_done_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_check_done_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_check_done();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_check_done_23_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_refine(const int *ic3, int io_reader);
 */
void CProxy_EnzoBlock::p_restart_refine(const int *ic3, int io_reader, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const int *ic3, int io_reader
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|io_reader;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|io_reader;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_restart_refine_marshall24(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_restart_refine_marshall24() {
  int epidx = CkRegisterEp("p_restart_refine(const int *ic3, int io_reader)",
      reinterpret_cast<CkCallFnPtr>(_call_p_restart_refine_marshall24), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_restart_refine_marshall24);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_restart_refine_marshall24);

  return epidx;
}

void CkIndex_EnzoBlock::_call_p_restart_refine_marshall24(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const int *ic3, int io_reader*/
  PUP::fromMem implP(impl_buf);
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  PUP::detail::TemporaryObjectHolder<int> io_reader;
  implP|io_reader;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  impl_obj->p_restart_refine(ic3, std::move(io_reader.t));
}
int CkIndex_EnzoBlock::_callmarshall_p_restart_refine_marshall24(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const int *ic3, int io_reader*/
  PUP::fromMem implP(impl_buf);
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  PUP::detail::TemporaryObjectHolder<int> io_reader;
  implP|io_reader;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  int *ic3=(int *)(impl_buf+impl_off_ic3);
  impl_obj->p_restart_refine(ic3, std::move(io_reader.t));
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_p_restart_refine_marshall24(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const int *ic3, int io_reader*/
  PUP::fromMem implP(impl_buf);
  int impl_off_ic3, impl_cnt_ic3;
  implP|impl_off_ic3;
  implP|impl_cnt_ic3;
  PUP::detail::TemporaryObjectHolder<int> io_reader;
  implP|io_reader;
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
  if (implDestP.hasComments()) implDestP.comment("io_reader");
  implDestP|io_reader;
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_restart_refine_24_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_set_data(EnzoMsgCheck* impl_msg);
 */
void CProxy_EnzoBlock::p_restart_set_data(EnzoMsgCheck* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_restart_set_data_EnzoMsgCheck(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_restart_set_data_EnzoMsgCheck() {
  int epidx = CkRegisterEp("p_restart_set_data(EnzoMsgCheck* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_restart_set_data_EnzoMsgCheck), CMessage_EnzoMsgCheck::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)EnzoMsgCheck::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_restart_set_data_EnzoMsgCheck(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_restart_set_data((EnzoMsgCheck*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_done();
 */
void CProxy_EnzoBlock::p_restart_done(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_restart_done_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_restart_done_void() {
  int epidx = CkRegisterEp("p_restart_done()",
      reinterpret_cast<CkCallFnPtr>(_call_p_restart_done_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_restart_done_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_restart_done();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_restart_done_26_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_accretion_end();
 */
void CProxy_EnzoBlock::p_method_accretion_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_method_accretion_end_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_method_accretion_end_void() {
  int epidx = CkRegisterEp("p_method_accretion_end()",
      reinterpret_cast<CkCallFnPtr>(_call_p_method_accretion_end_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_method_accretion_end_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_method_accretion_end();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_method_accretion_end_27_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_cg_matvec();
 */
void CProxy_EnzoBlock::p_solver_cg_matvec(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_cg_matvec_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_cg_matvec_void() {
  int epidx = CkRegisterEp("p_solver_cg_matvec()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_cg_matvec_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_cg_matvec_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_cg_matvec();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_cg_matvec_28_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_cg_loop_0a(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_0a_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_cg_loop_0a_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_cg_loop_0a(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_cg_loop_0a_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_cg_loop_0a_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_cg_loop_0a((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_cg_loop_0b(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_0b_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_cg_loop_0b_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_cg_loop_0b(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_cg_loop_0b_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_cg_loop_0b_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_cg_loop_0b((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_cg_shift_1(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_shift_1_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_cg_shift_1_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_cg_shift_1(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_cg_shift_1_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_cg_shift_1_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_cg_shift_1((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_cg_loop_2();
 */
void CProxy_EnzoBlock::p_solver_cg_loop_2(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_cg_loop_2_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_cg_loop_2_void() {
  int epidx = CkRegisterEp("p_solver_cg_loop_2()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_cg_loop_2_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_cg_loop_2_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_cg_loop_2();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_cg_loop_2_32_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_cg_loop_3(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_3_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_cg_loop_3_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_cg_loop_3(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_cg_loop_3_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_cg_loop_3_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_cg_loop_3((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_cg_loop_5(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_5_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_cg_loop_5_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_cg_loop_5(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_cg_loop_5_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_cg_loop_5_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_cg_loop_5((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_bicgstab_start_1(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_start_1_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_bicgstab_start_1_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_bicgstab_start_1(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_bicgstab_start_1_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_bicgstab_start_1_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_bicgstab_start_1((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_bicgstab_start_3(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_start_3_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_bicgstab_start_3_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_bicgstab_start_3(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_bicgstab_start_3_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_bicgstab_start_3_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_bicgstab_start_3((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_5_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_bicgstab_loop_5_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_bicgstab_loop_5_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_bicgstab_loop_5_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_bicgstab_loop_5((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_11_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_bicgstab_loop_11_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_bicgstab_loop_11_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_bicgstab_loop_11_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_bicgstab_loop_11((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_13_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_bicgstab_loop_13_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_bicgstab_loop_13_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_bicgstab_loop_13_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_bicgstab_loop_13((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_15_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_bicgstab_loop_15_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_bicgstab_loop_15_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_bicgstab_loop_15_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_bicgstab_loop_15((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_2();
 */
void CProxy_EnzoBlock::p_solver_bicgstab_loop_2(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_2_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_bicgstab_loop_2_void() {
  int epidx = CkRegisterEp("p_solver_bicgstab_loop_2()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_bicgstab_loop_2_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_bicgstab_loop_2_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_bicgstab_loop_2();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_bicgstab_loop_2_41_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_3();
 */
void CProxy_EnzoBlock::p_solver_bicgstab_loop_3(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_3_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_bicgstab_loop_3_void() {
  int epidx = CkRegisterEp("p_solver_bicgstab_loop_3()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_bicgstab_loop_3_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_bicgstab_loop_3_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_bicgstab_loop_3();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_bicgstab_loop_3_42_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_8();
 */
void CProxy_EnzoBlock::p_solver_bicgstab_loop_8(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_8_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_bicgstab_loop_8_void() {
  int epidx = CkRegisterEp("p_solver_bicgstab_loop_8()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_bicgstab_loop_8_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_bicgstab_loop_8_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_bicgstab_loop_8();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_bicgstab_loop_8_43_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_9();
 */
void CProxy_EnzoBlock::p_solver_bicgstab_loop_9(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_9_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_bicgstab_loop_9_void() {
  int epidx = CkRegisterEp("p_solver_bicgstab_loop_9()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_bicgstab_loop_9_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_bicgstab_loop_9_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_bicgstab_loop_9();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_bicgstab_loop_9_44_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
 */
void CProxy_EnzoBlock::p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_dot, impl_cnt_dot;
  impl_off_dot=impl_off=CK_ALIGN(impl_off,sizeof(long double));
  impl_off+=(impl_cnt_dot=sizeof(long double)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
    implP|iter;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
    implP|iter;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_dot,dot,impl_cnt_dot);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_dot_recv_parent_marshall45(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_dot_recv_parent_marshall45() {
  int epidx = CkRegisterEp("p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter)",
      reinterpret_cast<CkCallFnPtr>(_call_p_dot_recv_parent_marshall45), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_dot_recv_parent_marshall45);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_dot_recv_parent_marshall45);

  return epidx;
}

void CkIndex_EnzoBlock::_call_p_dot_recv_parent_marshall45(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_dot, impl_cnt_dot;
  implP|impl_off_dot;
  implP|impl_cnt_dot;
  PUP::detail::TemporaryObjectHolder<std::vector<int>> isa;
  implP|isa;
  PUP::detail::TemporaryObjectHolder<int> i_function;
  implP|i_function;
  PUP::detail::TemporaryObjectHolder<int> iter;
  implP|iter;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  long double *dot=(long double *)(impl_buf+impl_off_dot);
  impl_obj->p_dot_recv_parent(std::move(n.t), dot, std::move(isa.t), std::move(i_function.t), std::move(iter.t));
}
int CkIndex_EnzoBlock::_callmarshall_p_dot_recv_parent_marshall45(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_dot, impl_cnt_dot;
  implP|impl_off_dot;
  implP|impl_cnt_dot;
  PUP::detail::TemporaryObjectHolder<std::vector<int>> isa;
  implP|isa;
  PUP::detail::TemporaryObjectHolder<int> i_function;
  implP|i_function;
  PUP::detail::TemporaryObjectHolder<int> iter;
  implP|iter;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  long double *dot=(long double *)(impl_buf+impl_off_dot);
  impl_obj->p_dot_recv_parent(std::move(n.t), dot, std::move(isa.t), std::move(i_function.t), std::move(iter.t));
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_p_dot_recv_parent_marshall45(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_dot, impl_cnt_dot;
  implP|impl_off_dot;
  implP|impl_cnt_dot;
  PUP::detail::TemporaryObjectHolder<std::vector<int>> isa;
  implP|isa;
  PUP::detail::TemporaryObjectHolder<int> i_function;
  implP|i_function;
  PUP::detail::TemporaryObjectHolder<int> iter;
  implP|iter;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  long double *dot=(long double *)(impl_buf+impl_off_dot);
  if (implDestP.hasComments()) implDestP.comment("n");
  implDestP|n;
  if (implDestP.hasComments()) implDestP.comment("dot");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*dot))<impl_cnt_dot;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|dot[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("isa");
  implDestP|isa;
  if (implDestP.hasComments()) implDestP.comment("i_function");
  implDestP|i_function;
  if (implDestP.hasComments()) implDestP.comment("iter");
  implDestP|iter;
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_dot_recv_parent_45_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
 */
void CProxy_EnzoBlock::p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const long double *dot, const std::vector<int> &isa, int i_function
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_dot, impl_cnt_dot;
  impl_off_dot=impl_off=CK_ALIGN(impl_off,sizeof(long double));
  impl_off+=(impl_cnt_dot=sizeof(long double)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_dot,dot,impl_cnt_dot);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_dot_recv_children_marshall46(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_dot_recv_children_marshall46() {
  int epidx = CkRegisterEp("p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function)",
      reinterpret_cast<CkCallFnPtr>(_call_p_dot_recv_children_marshall46), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_dot_recv_children_marshall46);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_dot_recv_children_marshall46);

  return epidx;
}

void CkIndex_EnzoBlock::_call_p_dot_recv_children_marshall46(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const long double *dot, const std::vector<int> &isa, int i_function*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_dot, impl_cnt_dot;
  implP|impl_off_dot;
  implP|impl_cnt_dot;
  PUP::detail::TemporaryObjectHolder<std::vector<int>> isa;
  implP|isa;
  PUP::detail::TemporaryObjectHolder<int> i_function;
  implP|i_function;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  long double *dot=(long double *)(impl_buf+impl_off_dot);
  impl_obj->p_dot_recv_children(std::move(n.t), dot, std::move(isa.t), std::move(i_function.t));
}
int CkIndex_EnzoBlock::_callmarshall_p_dot_recv_children_marshall46(char* impl_buf, void* impl_obj_void) {
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int n, const long double *dot, const std::vector<int> &isa, int i_function*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_dot, impl_cnt_dot;
  implP|impl_off_dot;
  implP|impl_cnt_dot;
  PUP::detail::TemporaryObjectHolder<std::vector<int>> isa;
  implP|isa;
  PUP::detail::TemporaryObjectHolder<int> i_function;
  implP|i_function;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  long double *dot=(long double *)(impl_buf+impl_off_dot);
  impl_obj->p_dot_recv_children(std::move(n.t), dot, std::move(isa.t), std::move(i_function.t));
  return implP.size();
}
void CkIndex_EnzoBlock::_marshallmessagepup_p_dot_recv_children_marshall46(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const long double *dot, const std::vector<int> &isa, int i_function*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_dot, impl_cnt_dot;
  implP|impl_off_dot;
  implP|impl_cnt_dot;
  PUP::detail::TemporaryObjectHolder<std::vector<int>> isa;
  implP|isa;
  PUP::detail::TemporaryObjectHolder<int> i_function;
  implP|i_function;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  long double *dot=(long double *)(impl_buf+impl_off_dot);
  if (implDestP.hasComments()) implDestP.comment("n");
  implDestP|n;
  if (implDestP.hasComments()) implDestP.comment("dot");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*dot))<impl_cnt_dot;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|dot[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("isa");
  implDestP|isa;
  if (implDestP.hasComments()) implDestP.comment("i_function");
  implDestP|i_function;
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_dot_recv_children_46_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
 */
void CProxy_EnzoBlock::p_solver_dd_restrict_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_restrict_recv_FieldMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_dd_restrict_recv_FieldMsg() {
  int epidx = CkRegisterEp("p_solver_dd_restrict_recv(FieldMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_dd_restrict_recv_FieldMsg), CMessage_FieldMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)FieldMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_dd_restrict_recv_FieldMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_dd_restrict_recv((FieldMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
 */
void CProxy_EnzoBlock::p_solver_dd_prolong_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_prolong_recv_FieldMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_dd_prolong_recv_FieldMsg() {
  int epidx = CkRegisterEp("p_solver_dd_prolong_recv(FieldMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_dd_prolong_recv_FieldMsg), CMessage_FieldMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)FieldMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_dd_prolong_recv_FieldMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_dd_prolong_recv((FieldMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_solve_coarse();
 */
void CProxy_EnzoBlock::p_solver_dd_solve_coarse(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_solve_coarse_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_dd_solve_coarse_void() {
  int epidx = CkRegisterEp("p_solver_dd_solve_coarse()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_dd_solve_coarse_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_dd_solve_coarse_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_dd_solve_coarse();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_dd_solve_coarse_49_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_solve_domain();
 */
void CProxy_EnzoBlock::p_solver_dd_solve_domain(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_solve_domain_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_dd_solve_domain_void() {
  int epidx = CkRegisterEp("p_solver_dd_solve_domain()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_dd_solve_domain_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_dd_solve_domain_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_dd_solve_domain();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_dd_solve_domain_50_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_last_smooth();
 */
void CProxy_EnzoBlock::p_solver_dd_last_smooth(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_last_smooth_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_dd_last_smooth_void() {
  int epidx = CkRegisterEp("p_solver_dd_last_smooth()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_dd_last_smooth_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_dd_last_smooth_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_dd_last_smooth();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_dd_last_smooth_51_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_dd_barrier(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_dd_barrier(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_dd_barrier_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_dd_barrier_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_dd_barrier(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_dd_barrier_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_dd_barrier_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_dd_barrier((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_dd_end(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_dd_end(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_dd_end_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_dd_end_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_dd_end(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_dd_end_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_dd_end_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_dd_end((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_jacobi_continue();
 */
void CProxy_EnzoBlock::p_solver_jacobi_continue(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_jacobi_continue_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_jacobi_continue_void() {
  int epidx = CkRegisterEp("p_solver_jacobi_continue()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_jacobi_continue_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_jacobi_continue_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_jacobi_continue();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_jacobi_continue_54_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_restrict();
 */
void CProxy_EnzoBlock::p_solver_mg0_restrict(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_restrict_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_mg0_restrict_void() {
  int epidx = CkRegisterEp("p_solver_mg0_restrict()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_mg0_restrict_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_mg0_restrict_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_mg0_restrict();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_mg0_restrict_55_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_solve_coarse();
 */
void CProxy_EnzoBlock::p_solver_mg0_solve_coarse(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_solve_coarse_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_mg0_solve_coarse_void() {
  int epidx = CkRegisterEp("p_solver_mg0_solve_coarse()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_mg0_solve_coarse_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_mg0_solve_coarse_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_mg0_solve_coarse();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_mg0_solve_coarse_56_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_post_smooth();
 */
void CProxy_EnzoBlock::p_solver_mg0_post_smooth(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_post_smooth_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_mg0_post_smooth_void() {
  int epidx = CkRegisterEp("p_solver_mg0_post_smooth()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_mg0_post_smooth_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_mg0_post_smooth_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_mg0_post_smooth();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_mg0_post_smooth_57_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_last_smooth();
 */
void CProxy_EnzoBlock::p_solver_mg0_last_smooth(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_last_smooth_void(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_mg0_last_smooth_void() {
  int epidx = CkRegisterEp("p_solver_mg0_last_smooth()",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_mg0_last_smooth_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_mg0_last_smooth_void(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_mg0_last_smooth();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoBlock::p_solver_mg0_last_smooth_58_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_mg0_begin_solve(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_mg0_begin_solve_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_mg0_begin_solve_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_mg0_begin_solve(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_mg0_begin_solve_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_mg0_begin_solve_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_mg0_begin_solve((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
 */
void CProxy_EnzoBlock::r_solver_mg0_barrier(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_mg0_barrier_CkReductionMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_r_solver_mg0_barrier_CkReductionMsg() {
  int epidx = CkRegisterEp("r_solver_mg0_barrier(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_solver_mg0_barrier_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_r_solver_mg0_barrier_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->r_solver_mg0_barrier((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
 */
void CProxy_EnzoBlock::p_solver_mg0_prolong_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_prolong_recv_FieldMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_mg0_prolong_recv_FieldMsg() {
  int epidx = CkRegisterEp("p_solver_mg0_prolong_recv(FieldMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_mg0_prolong_recv_FieldMsg), CMessage_FieldMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)FieldMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_mg0_prolong_recv_FieldMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_mg0_prolong_recv((FieldMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
 */
void CProxy_EnzoBlock::p_solver_mg0_restrict_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_restrict_recv_FieldMsg(),0);
}

// Entry point registration function
int CkIndex_EnzoBlock::reg_p_solver_mg0_restrict_recv_FieldMsg() {
  int epidx = CkRegisterEp("p_solver_mg0_restrict_recv(FieldMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_solver_mg0_restrict_recv_FieldMsg), CMessage_FieldMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)FieldMsg::ckDebugPup);
  return epidx;
}

void CkIndex_EnzoBlock::_call_p_solver_mg0_restrict_recv_FieldMsg(void* impl_msg, void* impl_obj_void)
{
  EnzoBlock* impl_obj = static_cast<EnzoBlock*>(impl_obj_void);
  impl_obj->p_solver_mg0_restrict_recv((FieldMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoBlock(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_EnzoBlock::reg_EnzoBlock_CkMigrateMessage() {
  int epidx = CkRegisterEp("EnzoBlock(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_EnzoBlock_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoBlock::_call_EnzoBlock_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  call_migration_constructor<EnzoBlock> c = impl_obj_void;
  c((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
void CProxySection_EnzoBlock::p_set_msg_refine(MsgRefine* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_set_msg_refine_MsgRefine(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_msg_check(EnzoMsgCheck* impl_msg);
 */
void CProxySection_EnzoBlock::p_set_msg_check(EnzoMsgCheck* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_set_msg_check_EnzoMsgCheck(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoBlock();
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_feedback_starss_end();
 */
void CProxySection_EnzoBlock::p_method_feedback_starss_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_feedback_starss_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_m1_closure_solve_transport_eqn();
 */
void CProxySection_EnzoBlock::p_method_m1_closure_solve_transport_eqn(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_m1_closure_solve_transport_eqn_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_m1_closure_set_global_averages_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_method_turbulence_end(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_method_turbulence_end(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_method_turbulence_end_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_initial_hdf5_recv(MsgInitial* impl_msg);
 */
void CProxySection_EnzoBlock::p_initial_hdf5_recv(MsgInitial* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_initial_hdf5_recv_MsgInitial(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_balance_migrate();
 */
void CProxySection_EnzoBlock::p_method_balance_migrate(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_balance_migrate_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_balance_done();
 */
void CProxySection_EnzoBlock::p_method_balance_done(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_balance_done_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_gravity_continue();
 */
void CProxySection_EnzoBlock::p_method_gravity_continue(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_gravity_continue_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_gravity_end();
 */
void CProxySection_EnzoBlock::p_method_gravity_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_gravity_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
 */
void CProxySection_EnzoBlock::p_method_infer_merge_masks(int n, const char *mask, const int *ic3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const char *mask, const int *ic3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_mask, impl_cnt_mask;
  impl_off_mask=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_mask=sizeof(char)*(n));
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_mask;
    implP|impl_cnt_mask;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_mask;
    implP|impl_cnt_mask;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_mask,mask,impl_cnt_mask);
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_merge_masks_marshall14(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_count_arrays(int count);
 */
void CProxySection_EnzoBlock::p_method_infer_count_arrays(int count, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int count
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|count;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|count;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_count_arrays_marshall15(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_request_data(const int *il3);
 */
void CProxySection_EnzoBlock::p_method_infer_request_data(const int *il3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const int *il3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_il3, impl_cnt_il3;
  impl_off_il3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_il3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_il3;
    implP|impl_cnt_il3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_il3;
    implP|impl_cnt_il3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_il3,il3,impl_cnt_il3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_request_data_marshall16(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_update(int n, const char *buffer, const int *il3);
 */
void CProxySection_EnzoBlock::p_method_infer_update(int n, const char *buffer, const int *il3, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const char *buffer, const int *il3
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_buffer, impl_cnt_buffer;
  impl_off_buffer=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_buffer=sizeof(char)*(n));
  int impl_off_il3, impl_cnt_il3;
  impl_off_il3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_il3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    implP|impl_off_il3;
    implP|impl_cnt_il3;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    implP|impl_off_il3;
    implP|impl_cnt_il3;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_buffer,buffer,impl_cnt_buffer);
  memcpy(impl_buf+impl_off_il3,il3,impl_cnt_il3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_update_marshall17(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_infer_exit();
 */
void CProxySection_EnzoBlock::p_method_infer_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_infer_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::p_method_fbnet_update_mesh(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_fbnet_update_mesh_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_fbnet_exit();
 */
void CProxySection_EnzoBlock::p_method_fbnet_exit(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_fbnet_exit_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
 */
void CProxySection_EnzoBlock::p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int num_files, const std::string &ordering, const std::string &name_dir
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)name_dir;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)name_dir;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_check_write_first_marshall21(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_write_next(int num_files, const std::string &ordering);
 */
void CProxySection_EnzoBlock::p_check_write_next(int num_files, const std::string &ordering, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int num_files, const std::string &ordering
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_check_write_next_marshall22(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_check_done();
 */
void CProxySection_EnzoBlock::p_check_done(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_check_done_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_refine(const int *ic3, int io_reader);
 */
void CProxySection_EnzoBlock::p_restart_refine(const int *ic3, int io_reader, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const int *ic3, int io_reader
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_ic3, impl_cnt_ic3;
  impl_off_ic3=impl_off=CK_ALIGN(impl_off,sizeof(int));
  impl_off+=(impl_cnt_ic3=sizeof(int)*(3));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|io_reader;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_ic3;
    implP|impl_cnt_ic3;
    implP|io_reader;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_ic3,ic3,impl_cnt_ic3);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_restart_refine_marshall24(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_set_data(EnzoMsgCheck* impl_msg);
 */
void CProxySection_EnzoBlock::p_restart_set_data(EnzoMsgCheck* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_restart_set_data_EnzoMsgCheck(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_done();
 */
void CProxySection_EnzoBlock::p_restart_done(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_restart_done_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_method_accretion_end();
 */
void CProxySection_EnzoBlock::p_method_accretion_end(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_method_accretion_end_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_cg_matvec();
 */
void CProxySection_EnzoBlock::p_solver_cg_matvec(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_cg_matvec_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_cg_loop_0a(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_0a_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_cg_loop_0b(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_0b_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_cg_shift_1(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_shift_1_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_cg_loop_2();
 */
void CProxySection_EnzoBlock::p_solver_cg_loop_2(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_cg_loop_2_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_cg_loop_3(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_3_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_cg_loop_5(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_cg_loop_5_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_bicgstab_start_1(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_start_1_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_bicgstab_start_3(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_start_3_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_5_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_11_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_13_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_bicgstab_loop_15_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_2();
 */
void CProxySection_EnzoBlock::p_solver_bicgstab_loop_2(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_2_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_3();
 */
void CProxySection_EnzoBlock::p_solver_bicgstab_loop_3(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_3_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_8();
 */
void CProxySection_EnzoBlock::p_solver_bicgstab_loop_8(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_8_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_bicgstab_loop_9();
 */
void CProxySection_EnzoBlock::p_solver_bicgstab_loop_9(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_bicgstab_loop_9_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
 */
void CProxySection_EnzoBlock::p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_dot, impl_cnt_dot;
  impl_off_dot=impl_off=CK_ALIGN(impl_off,sizeof(long double));
  impl_off+=(impl_cnt_dot=sizeof(long double)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
    implP|iter;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
    implP|iter;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_dot,dot,impl_cnt_dot);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_dot_recv_parent_marshall45(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
 */
void CProxySection_EnzoBlock::p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int n, const long double *dot, const std::vector<int> &isa, int i_function
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_dot, impl_cnt_dot;
  impl_off_dot=impl_off=CK_ALIGN(impl_off,sizeof(long double));
  impl_off+=(impl_cnt_dot=sizeof(long double)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_dot;
    implP|impl_cnt_dot;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::vector<int>>::type>::type &)isa;
    implP|i_function;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_dot,dot,impl_cnt_dot);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_dot_recv_children_marshall46(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
 */
void CProxySection_EnzoBlock::p_solver_dd_restrict_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_restrict_recv_FieldMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
 */
void CProxySection_EnzoBlock::p_solver_dd_prolong_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_prolong_recv_FieldMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_solve_coarse();
 */
void CProxySection_EnzoBlock::p_solver_dd_solve_coarse(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_solve_coarse_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_solve_domain();
 */
void CProxySection_EnzoBlock::p_solver_dd_solve_domain(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_solve_domain_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_dd_last_smooth();
 */
void CProxySection_EnzoBlock::p_solver_dd_last_smooth(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_dd_last_smooth_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_dd_barrier(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_dd_barrier(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_dd_barrier_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_dd_end(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_dd_end(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_dd_end_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_jacobi_continue();
 */
void CProxySection_EnzoBlock::p_solver_jacobi_continue(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_jacobi_continue_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_restrict();
 */
void CProxySection_EnzoBlock::p_solver_mg0_restrict(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_restrict_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_solve_coarse();
 */
void CProxySection_EnzoBlock::p_solver_mg0_solve_coarse(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_solve_coarse_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_post_smooth();
 */
void CProxySection_EnzoBlock::p_solver_mg0_post_smooth(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_post_smooth_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_last_smooth();
 */
void CProxySection_EnzoBlock::p_solver_mg0_last_smooth(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_last_smooth_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_mg0_begin_solve(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_mg0_begin_solve_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
 */
void CProxySection_EnzoBlock::r_solver_mg0_barrier(CkReductionMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_r_solver_mg0_barrier_CkReductionMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
 */
void CProxySection_EnzoBlock::p_solver_mg0_prolong_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_prolong_recv_FieldMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
 */
void CProxySection_EnzoBlock::p_solver_mg0_restrict_recv(FieldMsg* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoBlock::idx_p_solver_mg0_restrict_recv_FieldMsg(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoBlock(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_EnzoBlock::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterArrayDimensions(__idx, -1);
  CkRegisterBase(__idx, CkIndex_Block::__idx);
  // REG: EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
  idx_EnzoBlock_marshall1();

  // REG: void p_set_msg_refine(MsgRefine* impl_msg);
  idx_p_set_msg_refine_MsgRefine();

  // REG: void p_set_msg_check(EnzoMsgCheck* impl_msg);
  idx_p_set_msg_check_EnzoMsgCheck();

  // REG: EnzoBlock();
  idx_EnzoBlock_void();
  CkRegisterDefaultCtor(__idx, idx_EnzoBlock_void());

  // REG: void p_method_feedback_starss_end();
  idx_p_method_feedback_starss_end_void();

  // REG: void p_method_m1_closure_solve_transport_eqn();
  idx_p_method_m1_closure_solve_transport_eqn_void();

  // REG: void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
  idx_p_method_m1_closure_set_global_averages_CkReductionMsg();

  // REG: void r_method_turbulence_end(CkReductionMsg* impl_msg);
  idx_r_method_turbulence_end_CkReductionMsg();

  // REG: void p_initial_hdf5_recv(MsgInitial* impl_msg);
  idx_p_initial_hdf5_recv_MsgInitial();

  // REG: void p_method_balance_migrate();
  idx_p_method_balance_migrate_void();

  // REG: void p_method_balance_done();
  idx_p_method_balance_done_void();

  // REG: void p_method_gravity_continue();
  idx_p_method_gravity_continue_void();

  // REG: void p_method_gravity_end();
  idx_p_method_gravity_end_void();

  // REG: void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
  idx_p_method_infer_merge_masks_marshall14();

  // REG: void p_method_infer_count_arrays(int count);
  idx_p_method_infer_count_arrays_marshall15();

  // REG: void p_method_infer_request_data(const int *il3);
  idx_p_method_infer_request_data_marshall16();

  // REG: void p_method_infer_update(int n, const char *buffer, const int *il3);
  idx_p_method_infer_update_marshall17();

  // REG: void p_method_infer_exit();
  idx_p_method_infer_exit_void();

  // REG: void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
  idx_p_method_fbnet_update_mesh_CkReductionMsg();

  // REG: void p_method_fbnet_exit();
  idx_p_method_fbnet_exit_void();

  // REG: void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
  idx_p_check_write_first_marshall21();

  // REG: void p_check_write_next(int num_files, const std::string &ordering);
  idx_p_check_write_next_marshall22();

  // REG: void p_check_done();
  idx_p_check_done_void();

  // REG: void p_restart_refine(const int *ic3, int io_reader);
  idx_p_restart_refine_marshall24();

  // REG: void p_restart_set_data(EnzoMsgCheck* impl_msg);
  idx_p_restart_set_data_EnzoMsgCheck();

  // REG: void p_restart_done();
  idx_p_restart_done_void();

  // REG: void p_method_accretion_end();
  idx_p_method_accretion_end_void();

  // REG: void p_solver_cg_matvec();
  idx_p_solver_cg_matvec_void();

  // REG: void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
  idx_r_solver_cg_loop_0a_CkReductionMsg();

  // REG: void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
  idx_r_solver_cg_loop_0b_CkReductionMsg();

  // REG: void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
  idx_r_solver_cg_shift_1_CkReductionMsg();

  // REG: void p_solver_cg_loop_2();
  idx_p_solver_cg_loop_2_void();

  // REG: void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
  idx_r_solver_cg_loop_3_CkReductionMsg();

  // REG: void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
  idx_r_solver_cg_loop_5_CkReductionMsg();

  // REG: void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
  idx_r_solver_bicgstab_start_1_CkReductionMsg();

  // REG: void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
  idx_r_solver_bicgstab_start_3_CkReductionMsg();

  // REG: void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
  idx_r_solver_bicgstab_loop_5_CkReductionMsg();

  // REG: void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
  idx_r_solver_bicgstab_loop_11_CkReductionMsg();

  // REG: void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
  idx_r_solver_bicgstab_loop_13_CkReductionMsg();

  // REG: void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
  idx_r_solver_bicgstab_loop_15_CkReductionMsg();

  // REG: void p_solver_bicgstab_loop_2();
  idx_p_solver_bicgstab_loop_2_void();

  // REG: void p_solver_bicgstab_loop_3();
  idx_p_solver_bicgstab_loop_3_void();

  // REG: void p_solver_bicgstab_loop_8();
  idx_p_solver_bicgstab_loop_8_void();

  // REG: void p_solver_bicgstab_loop_9();
  idx_p_solver_bicgstab_loop_9_void();

  // REG: void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
  idx_p_dot_recv_parent_marshall45();

  // REG: void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
  idx_p_dot_recv_children_marshall46();

  // REG: void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
  idx_p_solver_dd_restrict_recv_FieldMsg();

  // REG: void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
  idx_p_solver_dd_prolong_recv_FieldMsg();

  // REG: void p_solver_dd_solve_coarse();
  idx_p_solver_dd_solve_coarse_void();

  // REG: void p_solver_dd_solve_domain();
  idx_p_solver_dd_solve_domain_void();

  // REG: void p_solver_dd_last_smooth();
  idx_p_solver_dd_last_smooth_void();

  // REG: void r_solver_dd_barrier(CkReductionMsg* impl_msg);
  idx_r_solver_dd_barrier_CkReductionMsg();

  // REG: void r_solver_dd_end(CkReductionMsg* impl_msg);
  idx_r_solver_dd_end_CkReductionMsg();

  // REG: void p_solver_jacobi_continue();
  idx_p_solver_jacobi_continue_void();

  // REG: void p_solver_mg0_restrict();
  idx_p_solver_mg0_restrict_void();

  // REG: void p_solver_mg0_solve_coarse();
  idx_p_solver_mg0_solve_coarse_void();

  // REG: void p_solver_mg0_post_smooth();
  idx_p_solver_mg0_post_smooth_void();

  // REG: void p_solver_mg0_last_smooth();
  idx_p_solver_mg0_last_smooth_void();

  // REG: void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
  idx_r_solver_mg0_begin_solve_CkReductionMsg();

  // REG: void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
  idx_r_solver_mg0_barrier_CkReductionMsg();

  // REG: void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
  idx_p_solver_mg0_prolong_recv_FieldMsg();

  // REG: void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
  idx_p_solver_mg0_restrict_recv_FieldMsg();

  // REG: EnzoBlock(CkMigrateMessage* impl_msg);
  idx_EnzoBlock_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_EnzoBlock_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: array IoEnzoReader: IoReader{
IoEnzoReader();
void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
void p_create_level(int level);
void p_init_level(int level);
void p_block_created();
void p_block_ready();
IoEnzoReader(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_IoEnzoReader::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CProxySection_IoEnzoReader::contribute(CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, userData, fragSize);
}

void CProxySection_IoEnzoReader::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, userData, fragSize);
}

template <typename T>
void CProxySection_IoEnzoReader::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, userData, fragSize);
}

void CProxySection_IoEnzoReader::contribute(CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, cb, userData, fragSize);
}

void CProxySection_IoEnzoReader::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, cb, userData, fragSize);
}

template <typename T>
void CProxySection_IoEnzoReader::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, cb, userData, fragSize);
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoReader();
 */
void CProxyElement_IoEnzoReader::insert(int onPE, const CkEntryOptions *impl_e_opts)
{ 
   void *impl_msg = CkAllocSysMsg(impl_e_opts);
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_IoEnzoReader::idx_IoEnzoReader_void(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
 */
void CProxyElement_IoEnzoReader::p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const std::string &impl_noname_0, const std::string &impl_noname_1, int level
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_1;
    implP|level;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_1;
    implP|level;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_init_root_marshall2(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_create_level(int level);
 */
void CProxyElement_IoEnzoReader::p_create_level(int level, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int level
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|level;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|level;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_create_level_marshall3(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_init_level(int level);
 */
void CProxyElement_IoEnzoReader::p_init_level(int level, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int level
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|level;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|level;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_init_level_marshall4(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_block_created();
 */
void CProxyElement_IoEnzoReader::p_block_created(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_block_created_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_block_ready();
 */
void CProxyElement_IoEnzoReader::p_block_ready(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_block_ready_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoReader(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoReader();
 */
CkArrayID CProxy_IoEnzoReader::ckNew(const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoEnzoReader::idx_IoEnzoReader_void(), opts);
  return gId;
}
void CProxy_IoEnzoReader::ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoEnzoReader::idx_IoEnzoReader_void(), _ck_array_creation_cb, opts, impl_msg);
}
CkArrayID CProxy_IoEnzoReader::ckNew(const int s1, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoEnzoReader::idx_IoEnzoReader_void(), opts);
  return gId;
}
void CProxy_IoEnzoReader::ckNew(const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoEnzoReader::idx_IoEnzoReader_void(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_IoEnzoReader::reg_IoEnzoReader_void() {
  int epidx = CkRegisterEp("IoEnzoReader()",
      reinterpret_cast<CkCallFnPtr>(_call_IoEnzoReader_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoEnzoReader::_call_IoEnzoReader_void(void* impl_msg, void* impl_obj_void)
{
  IoEnzoReader* impl_obj = static_cast<IoEnzoReader*>(impl_obj_void);
  new (impl_obj_void) IoEnzoReader();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
 */
void CProxy_IoEnzoReader::p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const std::string &impl_noname_0, const std::string &impl_noname_1, int level
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_1;
    implP|level;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_1;
    implP|level;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_IoEnzoReader::idx_p_init_root_marshall2(),0);
}

// Entry point registration function
int CkIndex_IoEnzoReader::reg_p_init_root_marshall2() {
  int epidx = CkRegisterEp("p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level)",
      reinterpret_cast<CkCallFnPtr>(_call_p_init_root_marshall2), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_init_root_marshall2);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_init_root_marshall2);

  return epidx;
}

void CkIndex_IoEnzoReader::_call_p_init_root_marshall2(void* impl_msg, void* impl_obj_void)
{
  IoEnzoReader* impl_obj = static_cast<IoEnzoReader*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const std::string &impl_noname_0, const std::string &impl_noname_1, int level*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<std::string> impl_noname_0;
  implP|impl_noname_0;
  PUP::detail::TemporaryObjectHolder<std::string> impl_noname_1;
  implP|impl_noname_1;
  PUP::detail::TemporaryObjectHolder<int> level;
  implP|level;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_init_root(std::move(impl_noname_0.t), std::move(impl_noname_1.t), std::move(level.t));
}
int CkIndex_IoEnzoReader::_callmarshall_p_init_root_marshall2(char* impl_buf, void* impl_obj_void) {
  IoEnzoReader* impl_obj = static_cast<IoEnzoReader*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const std::string &impl_noname_0, const std::string &impl_noname_1, int level*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<std::string> impl_noname_0;
  implP|impl_noname_0;
  PUP::detail::TemporaryObjectHolder<std::string> impl_noname_1;
  implP|impl_noname_1;
  PUP::detail::TemporaryObjectHolder<int> level;
  implP|level;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_init_root(std::move(impl_noname_0.t), std::move(impl_noname_1.t), std::move(level.t));
  return implP.size();
}
void CkIndex_IoEnzoReader::_marshallmessagepup_p_init_root_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const std::string &impl_noname_0, const std::string &impl_noname_1, int level*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<std::string> impl_noname_0;
  implP|impl_noname_0;
  PUP::detail::TemporaryObjectHolder<std::string> impl_noname_1;
  implP|impl_noname_1;
  PUP::detail::TemporaryObjectHolder<int> level;
  implP|level;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_1");
  implDestP|impl_noname_1;
  if (implDestP.hasComments()) implDestP.comment("level");
  implDestP|level;
}
PUPable_def(SINGLE_ARG(Closure_IoEnzoReader::p_init_root_2_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_create_level(int level);
 */
void CProxy_IoEnzoReader::p_create_level(int level, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int level
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|level;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|level;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_IoEnzoReader::idx_p_create_level_marshall3(),0);
}

// Entry point registration function
int CkIndex_IoEnzoReader::reg_p_create_level_marshall3() {
  int epidx = CkRegisterEp("p_create_level(int level)",
      reinterpret_cast<CkCallFnPtr>(_call_p_create_level_marshall3), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_create_level_marshall3);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_create_level_marshall3);

  return epidx;
}

void CkIndex_IoEnzoReader::_call_p_create_level_marshall3(void* impl_msg, void* impl_obj_void)
{
  IoEnzoReader* impl_obj = static_cast<IoEnzoReader*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int level*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> level;
  implP|level;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_create_level(std::move(level.t));
}
int CkIndex_IoEnzoReader::_callmarshall_p_create_level_marshall3(char* impl_buf, void* impl_obj_void) {
  IoEnzoReader* impl_obj = static_cast<IoEnzoReader*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int level*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> level;
  implP|level;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_create_level(std::move(level.t));
  return implP.size();
}
void CkIndex_IoEnzoReader::_marshallmessagepup_p_create_level_marshall3(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int level*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> level;
  implP|level;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("level");
  implDestP|level;
}
PUPable_def(SINGLE_ARG(Closure_IoEnzoReader::p_create_level_3_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_init_level(int level);
 */
void CProxy_IoEnzoReader::p_init_level(int level, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int level
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|level;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|level;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_IoEnzoReader::idx_p_init_level_marshall4(),0);
}

// Entry point registration function
int CkIndex_IoEnzoReader::reg_p_init_level_marshall4() {
  int epidx = CkRegisterEp("p_init_level(int level)",
      reinterpret_cast<CkCallFnPtr>(_call_p_init_level_marshall4), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_init_level_marshall4);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_init_level_marshall4);

  return epidx;
}

void CkIndex_IoEnzoReader::_call_p_init_level_marshall4(void* impl_msg, void* impl_obj_void)
{
  IoEnzoReader* impl_obj = static_cast<IoEnzoReader*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int level*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> level;
  implP|level;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_init_level(std::move(level.t));
}
int CkIndex_IoEnzoReader::_callmarshall_p_init_level_marshall4(char* impl_buf, void* impl_obj_void) {
  IoEnzoReader* impl_obj = static_cast<IoEnzoReader*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int level*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> level;
  implP|level;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_init_level(std::move(level.t));
  return implP.size();
}
void CkIndex_IoEnzoReader::_marshallmessagepup_p_init_level_marshall4(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int level*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> level;
  implP|level;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("level");
  implDestP|level;
}
PUPable_def(SINGLE_ARG(Closure_IoEnzoReader::p_init_level_4_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_block_created();
 */
void CProxy_IoEnzoReader::p_block_created(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_IoEnzoReader::idx_p_block_created_void(),0);
}

// Entry point registration function
int CkIndex_IoEnzoReader::reg_p_block_created_void() {
  int epidx = CkRegisterEp("p_block_created()",
      reinterpret_cast<CkCallFnPtr>(_call_p_block_created_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoEnzoReader::_call_p_block_created_void(void* impl_msg, void* impl_obj_void)
{
  IoEnzoReader* impl_obj = static_cast<IoEnzoReader*>(impl_obj_void);
  impl_obj->p_block_created();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_IoEnzoReader::p_block_created_5_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_block_ready();
 */
void CProxy_IoEnzoReader::p_block_ready(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_IoEnzoReader::idx_p_block_ready_void(),0);
}

// Entry point registration function
int CkIndex_IoEnzoReader::reg_p_block_ready_void() {
  int epidx = CkRegisterEp("p_block_ready()",
      reinterpret_cast<CkCallFnPtr>(_call_p_block_ready_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoEnzoReader::_call_p_block_ready_void(void* impl_msg, void* impl_obj_void)
{
  IoEnzoReader* impl_obj = static_cast<IoEnzoReader*>(impl_obj_void);
  impl_obj->p_block_ready();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_IoEnzoReader::p_block_ready_6_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoReader(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_IoEnzoReader::reg_IoEnzoReader_CkMigrateMessage() {
  int epidx = CkRegisterEp("IoEnzoReader(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_IoEnzoReader_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoEnzoReader::_call_IoEnzoReader_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  call_migration_constructor<IoEnzoReader> c = impl_obj_void;
  c((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoReader();
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
 */
void CProxySection_IoEnzoReader::p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const std::string &impl_noname_0, const std::string &impl_noname_1, int level
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_1;
    implP|level;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_0;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_1;
    implP|level;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_init_root_marshall2(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_create_level(int level);
 */
void CProxySection_IoEnzoReader::p_create_level(int level, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int level
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|level;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|level;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_create_level_marshall3(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_init_level(int level);
 */
void CProxySection_IoEnzoReader::p_init_level(int level, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: int level
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|level;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|level;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_init_level_marshall4(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_block_created();
 */
void CProxySection_IoEnzoReader::p_block_created(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_block_created_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_block_ready();
 */
void CProxySection_IoEnzoReader::p_block_ready(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoReader::idx_p_block_ready_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoReader(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_IoEnzoReader::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterArrayDimensions(__idx, 1);
  CkRegisterBase(__idx, CkIndex_IoReader::__idx);
  // REG: IoEnzoReader();
  idx_IoEnzoReader_void();
  CkRegisterDefaultCtor(__idx, idx_IoEnzoReader_void());

  // REG: void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
  idx_p_init_root_marshall2();

  // REG: void p_create_level(int level);
  idx_p_create_level_marshall3();

  // REG: void p_init_level(int level);
  idx_p_init_level_marshall4();

  // REG: void p_block_created();
  idx_p_block_created_void();

  // REG: void p_block_ready();
  idx_p_block_ready_void();

  // REG: IoEnzoReader(CkMigrateMessage* impl_msg);
  idx_IoEnzoReader_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_IoEnzoReader_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: array IoEnzoWriter: IoWriter{
IoEnzoWriter();
IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
void p_write(EnzoMsgCheck* impl_msg);
IoEnzoWriter(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_IoEnzoWriter::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CProxySection_IoEnzoWriter::contribute(CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, userData, fragSize);
}

void CProxySection_IoEnzoWriter::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, userData, fragSize);
}

template <typename T>
void CProxySection_IoEnzoWriter::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, userData, fragSize);
}

void CProxySection_IoEnzoWriter::contribute(CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, cb, userData, fragSize);
}

void CProxySection_IoEnzoWriter::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, cb, userData, fragSize);
}

template <typename T>
void CProxySection_IoEnzoWriter::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, cb, userData, fragSize);
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoWriter();
 */
void CProxyElement_IoEnzoWriter::insert(int onPE, const CkEntryOptions *impl_e_opts)
{ 
   void *impl_msg = CkAllocSysMsg(impl_e_opts);
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_IoEnzoWriter::idx_IoEnzoWriter_void(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
 */
void CProxyElement_IoEnzoWriter::insert(int num_files, const std::string &ordering, int monitor_iter, int onPE, const CkEntryOptions *impl_e_opts)
{ 
   //Marshall: int num_files, const std::string &ordering, int monitor_iter
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
  }
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_IoEnzoWriter::idx_IoEnzoWriter_marshall2(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_write(EnzoMsgCheck* impl_msg);
 */
void CProxyElement_IoEnzoWriter::p_write(EnzoMsgCheck* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoWriter::idx_p_write_EnzoMsgCheck(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoWriter(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoWriter();
 */
CkArrayID CProxy_IoEnzoWriter::ckNew(const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoEnzoWriter::idx_IoEnzoWriter_void(), opts);
  return gId;
}
void CProxy_IoEnzoWriter::ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoEnzoWriter::idx_IoEnzoWriter_void(), _ck_array_creation_cb, opts, impl_msg);
}
CkArrayID CProxy_IoEnzoWriter::ckNew(const int s1, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoEnzoWriter::idx_IoEnzoWriter_void(), opts);
  return gId;
}
void CProxy_IoEnzoWriter::ckNew(const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoEnzoWriter::idx_IoEnzoWriter_void(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_IoEnzoWriter::reg_IoEnzoWriter_void() {
  int epidx = CkRegisterEp("IoEnzoWriter()",
      reinterpret_cast<CkCallFnPtr>(_call_IoEnzoWriter_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoEnzoWriter::_call_IoEnzoWriter_void(void* impl_msg, void* impl_obj_void)
{
  IoEnzoWriter* impl_obj = static_cast<IoEnzoWriter*>(impl_obj_void);
  new (impl_obj_void) IoEnzoWriter();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
 */
CkArrayID CProxy_IoEnzoWriter::ckNew(int num_files, const std::string &ordering, int monitor_iter, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int num_files, const std::string &ordering, int monitor_iter
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
  }
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoEnzoWriter::idx_IoEnzoWriter_marshall2(), opts);
  return gId;
}
void CProxy_IoEnzoWriter::ckNew(int num_files, const std::string &ordering, int monitor_iter, const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int num_files, const std::string &ordering, int monitor_iter
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
  }
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoEnzoWriter::idx_IoEnzoWriter_marshall2(), _ck_array_creation_cb, opts, impl_msg);
}
CkArrayID CProxy_IoEnzoWriter::ckNew(int num_files, const std::string &ordering, int monitor_iter, const int s1, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int num_files, const std::string &ordering, int monitor_iter
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
  }
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_IoEnzoWriter::idx_IoEnzoWriter_marshall2(), opts);
  return gId;
}
void CProxy_IoEnzoWriter::ckNew(int num_files, const std::string &ordering, int monitor_iter, const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int num_files, const std::string &ordering, int monitor_iter
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|num_files;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)ordering;
    implP|monitor_iter;
  }
  CkArrayOptions opts(s1);
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_IoEnzoWriter::idx_IoEnzoWriter_marshall2(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_IoEnzoWriter::reg_IoEnzoWriter_marshall2() {
  int epidx = CkRegisterEp("IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter)",
      reinterpret_cast<CkCallFnPtr>(_call_IoEnzoWriter_marshall2), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_IoEnzoWriter_marshall2);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_IoEnzoWriter_marshall2);

  return epidx;
}

void CkIndex_IoEnzoWriter::_call_IoEnzoWriter_marshall2(void* impl_msg, void* impl_obj_void)
{
  IoEnzoWriter* impl_obj = static_cast<IoEnzoWriter*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int num_files, const std::string &ordering, int monitor_iter*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> num_files;
  implP|num_files;
  PUP::detail::TemporaryObjectHolder<std::string> ordering;
  implP|ordering;
  PUP::detail::TemporaryObjectHolder<int> monitor_iter;
  implP|monitor_iter;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) IoEnzoWriter(std::move(num_files.t), std::move(ordering.t), std::move(monitor_iter.t));
}
int CkIndex_IoEnzoWriter::_callmarshall_IoEnzoWriter_marshall2(char* impl_buf, void* impl_obj_void) {
  IoEnzoWriter* impl_obj = static_cast<IoEnzoWriter*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int num_files, const std::string &ordering, int monitor_iter*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> num_files;
  implP|num_files;
  PUP::detail::TemporaryObjectHolder<std::string> ordering;
  implP|ordering;
  PUP::detail::TemporaryObjectHolder<int> monitor_iter;
  implP|monitor_iter;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) IoEnzoWriter(std::move(num_files.t), std::move(ordering.t), std::move(monitor_iter.t));
  return implP.size();
}
void CkIndex_IoEnzoWriter::_marshallmessagepup_IoEnzoWriter_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int num_files, const std::string &ordering, int monitor_iter*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> num_files;
  implP|num_files;
  PUP::detail::TemporaryObjectHolder<std::string> ordering;
  implP|ordering;
  PUP::detail::TemporaryObjectHolder<int> monitor_iter;
  implP|monitor_iter;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("num_files");
  implDestP|num_files;
  if (implDestP.hasComments()) implDestP.comment("ordering");
  implDestP|ordering;
  if (implDestP.hasComments()) implDestP.comment("monitor_iter");
  implDestP|monitor_iter;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_write(EnzoMsgCheck* impl_msg);
 */
void CProxy_IoEnzoWriter::p_write(EnzoMsgCheck* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_IoEnzoWriter::idx_p_write_EnzoMsgCheck(),0);
}

// Entry point registration function
int CkIndex_IoEnzoWriter::reg_p_write_EnzoMsgCheck() {
  int epidx = CkRegisterEp("p_write(EnzoMsgCheck* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_p_write_EnzoMsgCheck), CMessage_EnzoMsgCheck::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)EnzoMsgCheck::ckDebugPup);
  return epidx;
}

void CkIndex_IoEnzoWriter::_call_p_write_EnzoMsgCheck(void* impl_msg, void* impl_obj_void)
{
  IoEnzoWriter* impl_obj = static_cast<IoEnzoWriter*>(impl_obj_void);
  impl_obj->p_write((EnzoMsgCheck*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoWriter(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_IoEnzoWriter::reg_IoEnzoWriter_CkMigrateMessage() {
  int epidx = CkRegisterEp("IoEnzoWriter(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_IoEnzoWriter_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_IoEnzoWriter::_call_IoEnzoWriter_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  call_migration_constructor<IoEnzoWriter> c = impl_obj_void;
  c((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoWriter();
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_write(EnzoMsgCheck* impl_msg);
 */
void CProxySection_IoEnzoWriter::p_write(EnzoMsgCheck* impl_msg) 
{
  ckCheck();
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_IoEnzoWriter::idx_p_write_EnzoMsgCheck(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: IoEnzoWriter(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_IoEnzoWriter::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterArrayDimensions(__idx, 1);
  CkRegisterBase(__idx, CkIndex_IoWriter::__idx);
  // REG: IoEnzoWriter();
  idx_IoEnzoWriter_void();
  CkRegisterDefaultCtor(__idx, idx_IoEnzoWriter_void());

  // REG: IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
  idx_IoEnzoWriter_marshall2();

  // REG: void p_write(EnzoMsgCheck* impl_msg);
  idx_p_write_EnzoMsgCheck();

  // REG: IoEnzoWriter(CkMigrateMessage* impl_msg);
  idx_IoEnzoWriter_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_IoEnzoWriter_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: array EnzoLevelArray: ArrayElement{
EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
void p_request_data();
void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
void p_done(const Index &impl_noname_4);
EnzoLevelArray(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_EnzoLevelArray::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CProxySection_EnzoLevelArray::contribute(CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, userData, fragSize);
}

void CProxySection_EnzoLevelArray::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, userData, fragSize);
}

template <typename T>
void CProxySection_EnzoLevelArray::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, userData, fragSize);
}

void CProxySection_EnzoLevelArray::contribute(CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(sid, cb, userData, fragSize);
}

void CProxySection_EnzoLevelArray::contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(dataSize, data, type, sid, cb, userData, fragSize);
}

template <typename T>
void CProxySection_EnzoLevelArray::contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData, int fragSize)
{
   CkArray *ckarr = CProxy_CkArray(sid.get_aid()).ckLocalBranch();
   CkMulticastMgr *mCastGrp = CProxy_CkMulticastMgr(ckarr->getmCastMgr()).ckLocalBranch();
   mCastGrp->contribute(data, type, sid, cb, userData, fragSize);
}

#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
 */
void CProxyElement_EnzoLevelArray::insert(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz, int onPE, const CkEntryOptions *impl_e_opts)
{ 
   //Marshall: const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_2;
    implP|level_base;
    implP|level_array;
    implP|level_infer;
    implP|nax;
    implP|nay;
    implP|naz;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_2;
    implP|level_base;
    implP|level_array;
    implP|level_infer;
    implP|nax;
    implP|nay;
    implP|naz;
  }
   UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
   ckInsert((CkArrayMessage *)impl_msg,CkIndex_EnzoLevelArray::idx_EnzoLevelArray_marshall1(),onPE);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_request_data();
 */
void CProxyElement_EnzoLevelArray::p_request_data(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoLevelArray::idx_p_request_data_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
 */
void CProxyElement_EnzoLevelArray::p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const Index &impl_noname_3, int nf, const enzo_float *field_data
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_field_data, impl_cnt_field_data;
  impl_off_field_data=impl_off=CK_ALIGN(impl_off,sizeof(enzo_float));
  impl_off+=(impl_cnt_field_data=sizeof(enzo_float)*(nf));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_3;
    implP|nf;
    implP|impl_off_field_data;
    implP|impl_cnt_field_data;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_3;
    implP|nf;
    implP|impl_off_field_data;
    implP|impl_cnt_field_data;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_field_data,field_data,impl_cnt_field_data);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoLevelArray::idx_p_transfer_data_marshall3(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_done(const Index &impl_noname_4);
 */
void CProxyElement_EnzoLevelArray::p_done(const Index &impl_noname_4, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const Index &impl_noname_4
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_4;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_4;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoLevelArray::idx_p_done_marshall4(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoLevelArray(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
 */
CkArrayID CProxy_EnzoLevelArray::ckNew(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz, const CkArrayOptions &opts, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_2;
    implP|level_base;
    implP|level_array;
    implP|level_infer;
    implP|nax;
    implP|nay;
    implP|naz;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_2;
    implP|level_base;
    implP|level_array;
    implP|level_infer;
    implP|nax;
    implP|nay;
    implP|naz;
  }
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkArrayID gId = ckCreateArray((CkArrayMessage *)impl_msg, CkIndex_EnzoLevelArray::idx_EnzoLevelArray_marshall1(), opts);
  return gId;
}
void CProxy_EnzoLevelArray::ckNew(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz, const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_2;
    implP|level_base;
    implP|level_array;
    implP|level_infer;
    implP|nax;
    implP|nay;
    implP|naz;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)impl_noname_2;
    implP|level_base;
    implP|level_array;
    implP|level_infer;
    implP|nax;
    implP|nay;
    implP|naz;
  }
  UsrToEnv(impl_msg)->setMsgtype(ArrayEltInitMsg);
  CkSendAsyncCreateArray(CkIndex_EnzoLevelArray::idx_EnzoLevelArray_marshall1(), _ck_array_creation_cb, opts, impl_msg);
}

// Entry point registration function
int CkIndex_EnzoLevelArray::reg_EnzoLevelArray_marshall1() {
  int epidx = CkRegisterEp("EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz)",
      reinterpret_cast<CkCallFnPtr>(_call_EnzoLevelArray_marshall1), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_EnzoLevelArray_marshall1);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_EnzoLevelArray_marshall1);

  return epidx;
}

void CkIndex_EnzoLevelArray::_call_EnzoLevelArray_marshall1(void* impl_msg, void* impl_obj_void)
{
  EnzoLevelArray* impl_obj = static_cast<EnzoLevelArray*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<std::string> impl_noname_2;
  implP|impl_noname_2;
  PUP::detail::TemporaryObjectHolder<int> level_base;
  implP|level_base;
  PUP::detail::TemporaryObjectHolder<int> level_array;
  implP|level_array;
  PUP::detail::TemporaryObjectHolder<int> level_infer;
  implP|level_infer;
  PUP::detail::TemporaryObjectHolder<int> nax;
  implP|nax;
  PUP::detail::TemporaryObjectHolder<int> nay;
  implP|nay;
  PUP::detail::TemporaryObjectHolder<int> naz;
  implP|naz;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) EnzoLevelArray(std::move(impl_noname_2.t), std::move(level_base.t), std::move(level_array.t), std::move(level_infer.t), std::move(nax.t), std::move(nay.t), std::move(naz.t));
}
int CkIndex_EnzoLevelArray::_callmarshall_EnzoLevelArray_marshall1(char* impl_buf, void* impl_obj_void) {
  EnzoLevelArray* impl_obj = static_cast<EnzoLevelArray*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<std::string> impl_noname_2;
  implP|impl_noname_2;
  PUP::detail::TemporaryObjectHolder<int> level_base;
  implP|level_base;
  PUP::detail::TemporaryObjectHolder<int> level_array;
  implP|level_array;
  PUP::detail::TemporaryObjectHolder<int> level_infer;
  implP|level_infer;
  PUP::detail::TemporaryObjectHolder<int> nax;
  implP|nax;
  PUP::detail::TemporaryObjectHolder<int> nay;
  implP|nay;
  PUP::detail::TemporaryObjectHolder<int> naz;
  implP|naz;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) EnzoLevelArray(std::move(impl_noname_2.t), std::move(level_base.t), std::move(level_array.t), std::move(level_infer.t), std::move(nax.t), std::move(nay.t), std::move(naz.t));
  return implP.size();
}
void CkIndex_EnzoLevelArray::_marshallmessagepup_EnzoLevelArray_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<std::string> impl_noname_2;
  implP|impl_noname_2;
  PUP::detail::TemporaryObjectHolder<int> level_base;
  implP|level_base;
  PUP::detail::TemporaryObjectHolder<int> level_array;
  implP|level_array;
  PUP::detail::TemporaryObjectHolder<int> level_infer;
  implP|level_infer;
  PUP::detail::TemporaryObjectHolder<int> nax;
  implP|nax;
  PUP::detail::TemporaryObjectHolder<int> nay;
  implP|nay;
  PUP::detail::TemporaryObjectHolder<int> naz;
  implP|naz;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_2");
  implDestP|impl_noname_2;
  if (implDestP.hasComments()) implDestP.comment("level_base");
  implDestP|level_base;
  if (implDestP.hasComments()) implDestP.comment("level_array");
  implDestP|level_array;
  if (implDestP.hasComments()) implDestP.comment("level_infer");
  implDestP|level_infer;
  if (implDestP.hasComments()) implDestP.comment("nax");
  implDestP|nax;
  if (implDestP.hasComments()) implDestP.comment("nay");
  implDestP|nay;
  if (implDestP.hasComments()) implDestP.comment("naz");
  implDestP|naz;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_request_data();
 */
void CProxy_EnzoLevelArray::p_request_data(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoLevelArray::idx_p_request_data_void(),0);
}

// Entry point registration function
int CkIndex_EnzoLevelArray::reg_p_request_data_void() {
  int epidx = CkRegisterEp("p_request_data()",
      reinterpret_cast<CkCallFnPtr>(_call_p_request_data_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoLevelArray::_call_p_request_data_void(void* impl_msg, void* impl_obj_void)
{
  EnzoLevelArray* impl_obj = static_cast<EnzoLevelArray*>(impl_obj_void);
  impl_obj->p_request_data();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_EnzoLevelArray::p_request_data_2_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
 */
void CProxy_EnzoLevelArray::p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const Index &impl_noname_3, int nf, const enzo_float *field_data
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_field_data, impl_cnt_field_data;
  impl_off_field_data=impl_off=CK_ALIGN(impl_off,sizeof(enzo_float));
  impl_off+=(impl_cnt_field_data=sizeof(enzo_float)*(nf));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_3;
    implP|nf;
    implP|impl_off_field_data;
    implP|impl_cnt_field_data;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_3;
    implP|nf;
    implP|impl_off_field_data;
    implP|impl_cnt_field_data;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_field_data,field_data,impl_cnt_field_data);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoLevelArray::idx_p_transfer_data_marshall3(),0);
}

// Entry point registration function
int CkIndex_EnzoLevelArray::reg_p_transfer_data_marshall3() {
  int epidx = CkRegisterEp("p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data)",
      reinterpret_cast<CkCallFnPtr>(_call_p_transfer_data_marshall3), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_transfer_data_marshall3);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_transfer_data_marshall3);

  return epidx;
}

void CkIndex_EnzoLevelArray::_call_p_transfer_data_marshall3(void* impl_msg, void* impl_obj_void)
{
  EnzoLevelArray* impl_obj = static_cast<EnzoLevelArray*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &impl_noname_3, int nf, const enzo_float *field_data*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_3;
  implP|impl_noname_3;
  PUP::detail::TemporaryObjectHolder<int> nf;
  implP|nf;
  int impl_off_field_data, impl_cnt_field_data;
  implP|impl_off_field_data;
  implP|impl_cnt_field_data;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  enzo_float *field_data=(enzo_float *)(impl_buf+impl_off_field_data);
  impl_obj->p_transfer_data(std::move(impl_noname_3.t), std::move(nf.t), field_data);
}
int CkIndex_EnzoLevelArray::_callmarshall_p_transfer_data_marshall3(char* impl_buf, void* impl_obj_void) {
  EnzoLevelArray* impl_obj = static_cast<EnzoLevelArray*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const Index &impl_noname_3, int nf, const enzo_float *field_data*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_3;
  implP|impl_noname_3;
  PUP::detail::TemporaryObjectHolder<int> nf;
  implP|nf;
  int impl_off_field_data, impl_cnt_field_data;
  implP|impl_off_field_data;
  implP|impl_cnt_field_data;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  enzo_float *field_data=(enzo_float *)(impl_buf+impl_off_field_data);
  impl_obj->p_transfer_data(std::move(impl_noname_3.t), std::move(nf.t), field_data);
  return implP.size();
}
void CkIndex_EnzoLevelArray::_marshallmessagepup_p_transfer_data_marshall3(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &impl_noname_3, int nf, const enzo_float *field_data*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_3;
  implP|impl_noname_3;
  PUP::detail::TemporaryObjectHolder<int> nf;
  implP|nf;
  int impl_off_field_data, impl_cnt_field_data;
  implP|impl_off_field_data;
  implP|impl_cnt_field_data;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  enzo_float *field_data=(enzo_float *)(impl_buf+impl_off_field_data);
  if (implDestP.hasComments()) implDestP.comment("impl_noname_3");
  implDestP|impl_noname_3;
  if (implDestP.hasComments()) implDestP.comment("nf");
  implDestP|nf;
  if (implDestP.hasComments()) implDestP.comment("field_data");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*field_data))<impl_cnt_field_data;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|field_data[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_EnzoLevelArray::p_transfer_data_3_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_done(const Index &impl_noname_4);
 */
void CProxy_EnzoLevelArray::p_done(const Index &impl_noname_4, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const Index &impl_noname_4
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_4;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_4;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckBroadcast(impl_amsg, CkIndex_EnzoLevelArray::idx_p_done_marshall4(),0);
}

// Entry point registration function
int CkIndex_EnzoLevelArray::reg_p_done_marshall4() {
  int epidx = CkRegisterEp("p_done(const Index &impl_noname_4)",
      reinterpret_cast<CkCallFnPtr>(_call_p_done_marshall4), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_done_marshall4);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_done_marshall4);

  return epidx;
}

void CkIndex_EnzoLevelArray::_call_p_done_marshall4(void* impl_msg, void* impl_obj_void)
{
  EnzoLevelArray* impl_obj = static_cast<EnzoLevelArray*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &impl_noname_4*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_4;
  implP|impl_noname_4;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_done(std::move(impl_noname_4.t));
}
int CkIndex_EnzoLevelArray::_callmarshall_p_done_marshall4(char* impl_buf, void* impl_obj_void) {
  EnzoLevelArray* impl_obj = static_cast<EnzoLevelArray*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const Index &impl_noname_4*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_4;
  implP|impl_noname_4;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_done(std::move(impl_noname_4.t));
  return implP.size();
}
void CkIndex_EnzoLevelArray::_marshallmessagepup_p_done_marshall4(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &impl_noname_4*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_4;
  implP|impl_noname_4;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_4");
  implDestP|impl_noname_4;
}
PUPable_def(SINGLE_ARG(Closure_EnzoLevelArray::p_done_4_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoLevelArray(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_EnzoLevelArray::reg_EnzoLevelArray_CkMigrateMessage() {
  int epidx = CkRegisterEp("EnzoLevelArray(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_EnzoLevelArray_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_EnzoLevelArray::_call_EnzoLevelArray_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  call_migration_constructor<EnzoLevelArray> c = impl_obj_void;
  c((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_request_data();
 */
void CProxySection_EnzoLevelArray::p_request_data(const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoLevelArray::idx_p_request_data_void(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
 */
void CProxySection_EnzoLevelArray::p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const Index &impl_noname_3, int nf, const enzo_float *field_data
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_field_data, impl_cnt_field_data;
  impl_off_field_data=impl_off=CK_ALIGN(impl_off,sizeof(enzo_float));
  impl_off+=(impl_cnt_field_data=sizeof(enzo_float)*(nf));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_3;
    implP|nf;
    implP|impl_off_field_data;
    implP|impl_cnt_field_data;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_3;
    implP|nf;
    implP|impl_off_field_data;
    implP|impl_cnt_field_data;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_field_data,field_data,impl_cnt_field_data);
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoLevelArray::idx_p_transfer_data_marshall3(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_done(const Index &impl_noname_4);
 */
void CProxySection_EnzoLevelArray::p_done(const Index &impl_noname_4, const CkEntryOptions *impl_e_opts) 
{
  ckCheck();
  //Marshall: const Index &impl_noname_4
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_4;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_4;
  }
  UsrToEnv(impl_msg)->setMsgtype(ForArrayEltMsg);
  CkArrayMessage *impl_amsg=(CkArrayMessage *)impl_msg;
  impl_amsg->array_setIfNotThere(CkArray_IfNotThere_buffer);
  ckSend(impl_amsg, CkIndex_EnzoLevelArray::idx_p_done_marshall4(),0);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: EnzoLevelArray(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_EnzoLevelArray::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeArray);
  CkRegisterArrayDimensions(__idx, -1);
  CkRegisterBase(__idx, CkIndex_ArrayElement::__idx);
  // REG: EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
  idx_EnzoLevelArray_marshall1();

  // REG: void p_request_data();
  idx_p_request_data_void();

  // REG: void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
  idx_p_transfer_data_marshall3();

  // REG: void p_done(const Index &impl_noname_4);
  idx_p_done_marshall4();

  // REG: EnzoLevelArray(CkMigrateMessage* impl_msg);
  idx_EnzoLevelArray_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_EnzoLevelArray_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
void _registerenzo(void)
{
  static int _done = 0; if(_done) return; _done = 1;
  _registerInitCall(register_method_turbulence,1);

  _registerInitCall(mutex_init,1);

  _registerInitCall(mutex_init_bcg_iter,1);

  CkRegisterReadonly("EnzoMsgCheck::counter","int",sizeof(EnzoMsgCheck::counter),(void *) &EnzoMsgCheck::counter,__xlater_roPup_EnzoMsgCheck_QColon__QColon_counter);

  CkRegisterReadonly("g_enzo_config","EnzoConfig",sizeof(g_enzo_config),(void *) &g_enzo_config,__xlater_roPup_g_enzo_config);

  CkRegisterReadonly("EnzoBlock::UseMinimumPressureSupport","int",sizeof(EnzoBlock::UseMinimumPressureSupport),(void *) &EnzoBlock::UseMinimumPressureSupport,__xlater_roPup_EnzoBlock_QColon__QColon_UseMinimumPressureSupport);

  CkRegisterReadonly("EnzoBlock::MinimumPressureSupportParameter","enzo_float",sizeof(EnzoBlock::MinimumPressureSupportParameter),(void *) &EnzoBlock::MinimumPressureSupportParameter,__xlater_roPup_EnzoBlock_QColon__QColon_MinimumPressureSupportParameter);

  CkRegisterReadonly("EnzoBlock::MultiSpecies","int",sizeof(EnzoBlock::MultiSpecies),(void *) &EnzoBlock::MultiSpecies,__xlater_roPup_EnzoBlock_QColon__QColon_MultiSpecies);

  CkRegisterReadonly("EnzoBlock::PressureFree","int",sizeof(EnzoBlock::PressureFree),(void *) &EnzoBlock::PressureFree,__xlater_roPup_EnzoBlock_QColon__QColon_PressureFree);

  CkRegisterReadonly("EnzoBlock::GravitationalConstant","enzo_float",sizeof(EnzoBlock::GravitationalConstant),(void *) &EnzoBlock::GravitationalConstant,__xlater_roPup_EnzoBlock_QColon__QColon_GravitationalConstant);

  CkRegisterReadonly("EnzoBlock::ProblemType","int",sizeof(EnzoBlock::ProblemType),(void *) &EnzoBlock::ProblemType,__xlater_roPup_EnzoBlock_QColon__QColon_ProblemType);

  CkRegisterReadonly("EnzoBlock::PPMFlatteningParameter","int",sizeof(EnzoBlock::PPMFlatteningParameter),(void *) &EnzoBlock::PPMFlatteningParameter,__xlater_roPup_EnzoBlock_QColon__QColon_PPMFlatteningParameter);

  CkRegisterReadonly("EnzoBlock::PPMDiffusionParameter","int",sizeof(EnzoBlock::PPMDiffusionParameter),(void *) &EnzoBlock::PPMDiffusionParameter,__xlater_roPup_EnzoBlock_QColon__QColon_PPMDiffusionParameter);

  CkRegisterReadonly("EnzoBlock::PPMSteepeningParameter","int",sizeof(EnzoBlock::PPMSteepeningParameter),(void *) &EnzoBlock::PPMSteepeningParameter,__xlater_roPup_EnzoBlock_QColon__QColon_PPMSteepeningParameter);

  CkRegisterReadonly("EnzoBlock::InitialRedshift","enzo_float",sizeof(EnzoBlock::InitialRedshift),(void *) &EnzoBlock::InitialRedshift,__xlater_roPup_EnzoBlock_QColon__QColon_InitialRedshift);

  CkRegisterReadonly("EnzoBlock::InitialTimeInCodeUnits","enzo_float",sizeof(EnzoBlock::InitialTimeInCodeUnits),(void *) &EnzoBlock::InitialTimeInCodeUnits,__xlater_roPup_EnzoBlock_QColon__QColon_InitialTimeInCodeUnits);

  CkRegisterReadonly("EnzoBlock::DomainLeftEdge","enzo_float",sizeof(EnzoBlock::DomainLeftEdge),(void *) &EnzoBlock::DomainLeftEdge,__xlater_roPup_EnzoBlock_QColon__QColon_DomainLeftEdge);

  CkRegisterReadonly("EnzoBlock::DomainRightEdge","enzo_float",sizeof(EnzoBlock::DomainRightEdge),(void *) &EnzoBlock::DomainRightEdge,__xlater_roPup_EnzoBlock_QColon__QColon_DomainRightEdge);

  CkRegisterReadonly("EnzoBlock::GridRank","int",sizeof(EnzoBlock::GridRank),(void *) &EnzoBlock::GridRank,__xlater_roPup_EnzoBlock_QColon__QColon_GridRank);

  CkRegisterReadonly("EnzoBlock::ghost_depth","int",sizeof(EnzoBlock::ghost_depth),(void *) &EnzoBlock::ghost_depth,__xlater_roPup_EnzoBlock_QColon__QColon_ghost_depth);

  CkRegisterReadonly("EnzoBlock::NumberOfBaryonFields","int",sizeof(EnzoBlock::NumberOfBaryonFields),(void *) &EnzoBlock::NumberOfBaryonFields,__xlater_roPup_EnzoBlock_QColon__QColon_NumberOfBaryonFields);

      PUPable_reg(EnzoBoundary);

      PUPable_reg(EnzoConfig);

      PUPable_reg(EnzoFactory);

      PUPable_reg(EnzoInitialGrackleTest);

      PUPable_reg(EnzoInitialAccretionTest);

      PUPable_reg(EnzoInitialBBTest);

      PUPable_reg(EnzoInitialBCenter);

      PUPable_reg(EnzoInitialBurkertBodenheimer);

      PUPable_reg(EnzoInitialCloud);

      PUPable_reg(EnzoInitialCollapse);

      PUPable_reg(EnzoInitialCosmology);

      PUPable_reg(EnzoInitialFeedbackTest);

      PUPable_reg(EnzoInitialHdf5);

      PUPable_reg(EnzoInitialImplosion2);

      PUPable_reg(EnzoInitialInclinedWave);

      PUPable_reg(EnzoInitialIsolatedGalaxy);

      PUPable_reg(EnzoInitialMergeSinksTest);

      PUPable_reg(EnzoInitialMusic);

      PUPable_reg(EnzoInitialPm);

      PUPable_reg(EnzoInitialPpmlTest);

      PUPable_reg(EnzoInitialSedovArray2);

      PUPable_reg(EnzoInitialSedovArray3);

      PUPable_reg(EnzoInitialSedovRandom);

      PUPable_reg(EnzoInitialShockTube);

      PUPable_reg(EnzoInitialShuCollapse);

      PUPable_reg(EnzoInitialSoup);

      PUPable_reg(EnzoInitialTurbulence);

      PUPable_reg(IoEnzoBlock);

      PUPable_reg(EnzoRefineMass);

      PUPable_reg(EnzoRefineParticleMass);

      PUPable_reg(EnzoRefineShock);

      PUPable_reg(EnzoComputeCoolingTime);

      PUPable_reg(EnzoComputePressure);

      PUPable_reg(EnzoComputeTemperature);

      PUPable_reg(EnzoComputeAcceleration);

      PUPable_reg(EnzoComputeCicInterp);

      PUPable_reg(EnzoMatrixLaplace);

      PUPable_reg(EnzoMatrixDiagonal);

      PUPable_reg(EnzoMatrixIdentity);

      PUPable_reg(EnzoMethodGrackle);

      PUPable_reg(EnzoMethodAccretion);

      PUPable_reg(EnzoMethodBackgroundAcceleration);

      PUPable_reg(EnzoMethodBondiHoyleAccretion);

      PUPable_reg(EnzoMethodCheck);

      PUPable_reg(EnzoMethodComovingExpansion);

      PUPable_reg(EnzoMethodCosmology);

      PUPable_reg(EnzoMethodDistributedFeedback);

      PUPable_reg(EnzoMethodFeedback);

      PUPable_reg(EnzoMethodFeedbackSTARSS);

      PUPable_reg(EnzoMethodM1Closure);

      PUPable_reg(EnzoMethodFluxAccretion);

      PUPable_reg(EnzoMethodGravity);

      PUPable_reg(EnzoMethodHeat);

      PUPable_reg(EnzoMethodInference);

      PUPable_reg(EnzoMethodFBNetDeposit);

      PUPable_reg(EnzoMethodMergeSinks);

      PUPable_reg(EnzoMethodMHDVlct);

      PUPable_reg(EnzoMethodPmDeposit);

      PUPable_reg(EnzoMethodPmUpdate);

      PUPable_reg(EnzoMethodPpm);

      PUPable_reg(EnzoMethodPpml);

      PUPable_reg(EnzoMethodBalance);

      PUPable_reg(EnzoMethodSinkMaker);

      PUPable_reg(EnzoMethodThresholdAccretion);

      PUPable_reg(EnzoMethodTurbulence);

      PUPable_reg(EnzoIntegrationQuanUpdate);

      PUPable_reg(EnzoReconstructorNN);

      #define _CHARMXI_CLASS_NAME EnzoReconstructorPLM<PLM_EnzoRKLimiter>
      PUPable_reg2(_CHARMXI_CLASS_NAME,"EnzoReconstructorPLM<PLM_EnzoRKLimiter>");
      #undef _CHARMXI_CLASS_NAME

      #define _CHARMXI_CLASS_NAME EnzoReconstructorPLM<PLM_AthenaLimiter>
      PUPable_reg2(_CHARMXI_CLASS_NAME,"EnzoReconstructorPLM<PLM_AthenaLimiter>");
      #undef _CHARMXI_CLASS_NAME

      PUPable_reg(EnzoMethodStarMaker);

      PUPable_reg(EnzoMethodStarMakerStochasticSF);

      PUPable_reg(EnzoMethodStarMakerSTARSS);

      PUPable_reg(EnzoBfieldMethodCT);

      PUPable_reg(EnzoObjectFeedbackSphere);

      PUPable_reg(EnzoPhysicsCosmology);

      PUPable_reg(EnzoPhysicsFluidProps);

      PUPable_reg(EnzoProblem);

      PUPable_reg(EnzoProlong);

      PUPable_reg(EnzoProlongMC1);

      PUPable_reg(EnzoProlongPoisson);

      PUPable_reg(EnzoRestrict);

      PUPable_reg(EnzoSolverCg);

      PUPable_reg(EnzoSolverDd);

      PUPable_reg(EnzoSolverDiagonal);

      PUPable_reg(EnzoSolverBiCgStab);

      PUPable_reg(EnzoSolverMg0);

      PUPable_reg(EnzoSolverJacobi);

      PUPable_reg(EnzoStopping);

      PUPable_reg(EnzoUnits);

/* REG: message EnzoMsgCheck;
*/
CMessage_EnzoMsgCheck::__register("EnzoMsgCheck", sizeof(EnzoMsgCheck),(CkPackFnPtr) EnzoMsgCheck::pack,(CkUnpackFnPtr) EnzoMsgCheck::unpack);

  _registermesh();

  CkRegisterReadonly("proxy_enzo_simulation","CProxy_EnzoSimulation",sizeof(proxy_enzo_simulation),(void *) &proxy_enzo_simulation,__xlater_roPup_proxy_enzo_simulation);

  CkRegisterReadonly("proxy_io_enzo_reader","CProxy_IoEnzoReader",sizeof(proxy_io_enzo_reader),(void *) &proxy_io_enzo_reader,__xlater_roPup_proxy_io_enzo_reader);

  CkRegisterReadonly("proxy_io_enzo_writer","CProxy_IoEnzoWriter",sizeof(proxy_io_enzo_writer),(void *) &proxy_io_enzo_writer,__xlater_roPup_proxy_io_enzo_writer);

  CkRegisterReadonly("proxy_level_array","CProxy_EnzoLevelArray",sizeof(proxy_level_array),(void *) &proxy_level_array,__xlater_roPup_proxy_level_array);

/* REG: group EnzoSimulation: Simulation{
EnzoSimulation(const char *filename, int n);
void r_startup_begun(CkReductionMsg* impl_msg);
void p_get_msg_refine(const Index &index);
void p_get_msg_check(const Index &index);
void r_method_balance_count(CkReductionMsg* impl_msg);
void p_method_balance_check();
void r_method_check_enter(CkReductionMsg* impl_msg);
void p_check_done();
void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
void p_infer_set_array_count(int count);
void p_infer_array_created();
void p_infer_done();
void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
void p_io_reader_created();
void p_restart_next_level();
void p_restart_level_created();
void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
void p_fbnet_concatenate_sphere_lists();
void p_fbnet_done();
EnzoSimulation(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_EnzoSimulation::__register("EnzoSimulation", sizeof(EnzoSimulation));

/* REG: array EnzoBlock: Block{
EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
void p_set_msg_refine(MsgRefine* impl_msg);
void p_set_msg_check(EnzoMsgCheck* impl_msg);
EnzoBlock();
void p_method_feedback_starss_end();
void p_method_m1_closure_solve_transport_eqn();
void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
void r_method_turbulence_end(CkReductionMsg* impl_msg);
void p_initial_hdf5_recv(MsgInitial* impl_msg);
void p_method_balance_migrate();
void p_method_balance_done();
void p_method_gravity_continue();
void p_method_gravity_end();
void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
void p_method_infer_count_arrays(int count);
void p_method_infer_request_data(const int *il3);
void p_method_infer_update(int n, const char *buffer, const int *il3);
void p_method_infer_exit();
void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
void p_method_fbnet_exit();
void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
void p_check_write_next(int num_files, const std::string &ordering);
void p_check_done();
void p_restart_refine(const int *ic3, int io_reader);
void p_restart_set_data(EnzoMsgCheck* impl_msg);
void p_restart_done();
void p_method_accretion_end();
void p_solver_cg_matvec();
void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
void p_solver_cg_loop_2();
void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
void p_solver_bicgstab_loop_2();
void p_solver_bicgstab_loop_3();
void p_solver_bicgstab_loop_8();
void p_solver_bicgstab_loop_9();
void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
void p_solver_dd_solve_coarse();
void p_solver_dd_solve_domain();
void p_solver_dd_last_smooth();
void r_solver_dd_barrier(CkReductionMsg* impl_msg);
void r_solver_dd_end(CkReductionMsg* impl_msg);
void p_solver_jacobi_continue();
void p_solver_mg0_restrict();
void p_solver_mg0_solve_coarse();
void p_solver_mg0_post_smooth();
void p_solver_mg0_last_smooth();
void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
EnzoBlock(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_EnzoBlock::__register("EnzoBlock", sizeof(EnzoBlock));

/* REG: array IoEnzoReader: IoReader{
IoEnzoReader();
void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
void p_create_level(int level);
void p_init_level(int level);
void p_block_created();
void p_block_ready();
IoEnzoReader(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_IoEnzoReader::__register("IoEnzoReader", sizeof(IoEnzoReader));

/* REG: array IoEnzoWriter: IoWriter{
IoEnzoWriter();
IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
void p_write(EnzoMsgCheck* impl_msg);
IoEnzoWriter(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_IoEnzoWriter::__register("IoEnzoWriter", sizeof(IoEnzoWriter));

/* REG: array EnzoLevelArray: ArrayElement{
EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
void p_request_data();
void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
void p_done(const Index &impl_noname_4);
EnzoLevelArray(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_EnzoLevelArray::__register("EnzoLevelArray", sizeof(EnzoLevelArray));

}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_EnzoSimulation::virtual_pup(PUP::er &p) {
    recursive_pup<EnzoSimulation>(dynamic_cast<EnzoSimulation*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_EnzoBlock::virtual_pup(PUP::er &p) {
    recursive_pup<EnzoBlock>(dynamic_cast<EnzoBlock*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_IoEnzoReader::virtual_pup(PUP::er &p) {
    recursive_pup<IoEnzoReader>(dynamic_cast<IoEnzoReader*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_IoEnzoWriter::virtual_pup(PUP::er &p) {
    recursive_pup<IoEnzoWriter>(dynamic_cast<IoEnzoWriter*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_EnzoLevelArray::virtual_pup(PUP::er &p) {
    recursive_pup<EnzoLevelArray>(dynamic_cast<EnzoLevelArray*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
