/* Generated code for Python module '__main__'
 * created by Nuitka version 0.7.7
 *
 * This code is in part copyright 2021 Kay Hayen.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "nuitka/prelude.h"

#include "nuitka/unfreezing.h"

#include "__helpers.h"

/* The "module___main__" is a Python object pointer of module type.
 *
 * Note: For full compatibility with CPython, every module variable access
 * needs to go through it except for cases where the module cannot possibly
 * have changed in the mean time.
 */

PyObject *module___main__;
PyDictObject *moduledict___main__;

/* The declarations of module constants used, if any. */
static PyObject *mod_consts[109];
#ifndef __NUITKA_NO_ASSERT__
static Py_hash_t mod_consts_hash[109];
#endif

static PyObject *module_filename_obj = NULL;

/* Indicator if this modules private constants were created yet. */
static bool constants_created = false;

/* Function to create module private constants. */
static void createModuleConstants(void) {
    if (constants_created == false) {
        loadConstantsBlob(&mod_consts[0], UNTRANSLATE("__main__"));
        constants_created = true;

#ifndef __NUITKA_NO_ASSERT__
        for(int i = 0; i < 109; i++) {
            mod_consts_hash[i] = DEEP_HASH(mod_consts[i]);
        }
#endif
    }
}

// We want to be able to initialize the "__main__" constants in any case.
#if 1
void createMainModuleConstants(void) {
    createModuleConstants();
}
#endif

/* Function to verify module private constants for non-corruption. */
#ifndef __NUITKA_NO_ASSERT__
void checkModuleConstants___main__(void) {
    // The module may not have been used at all, then ignore this.
    if (constants_created == false) return;

    for(int i = 0; i < 109; i++) {
        assert(mod_consts_hash[i] == DEEP_HASH(mod_consts[i]));
        CHECK_OBJECT_DEEP(mod_consts[i]);
    }
}
#endif

// The module code objects.
static PyCodeObject *codeobj_ccd1fb807920dfdbbfd5c88e344983c6;
static PyCodeObject *codeobj_a4a3d272fc8190d3932415e3dda5345c;
static PyCodeObject *codeobj_7fdb2c580e4cfeb4386220f78b54ec88;
/* For use in "MainProgram.c". */
PyCodeObject *codeobj_main = NULL;

static void createModuleCodeObjects(void) {
    module_filename_obj = mod_consts[88]; CHECK_OBJECT(module_filename_obj);
    codeobj_ccd1fb807920dfdbbfd5c88e344983c6 = MAKE_CODEOBJECT(module_filename_obj, 71, CO_OPTIMIZED | CO_NEWLOCALS | CO_NOFREE, mod_consts[105], mod_consts[106], NULL, 1, 0, 0);
    codeobj_a4a3d272fc8190d3932415e3dda5345c = MAKE_CODEOBJECT(module_filename_obj, 1, CO_NOFREE, mod_consts[107], NULL, NULL, 0, 0, 0);
    codeobj_7fdb2c580e4cfeb4386220f78b54ec88 = MAKE_CODEOBJECT(module_filename_obj, 23, CO_OPTIMIZED | CO_NEWLOCALS | CO_NOFREE, mod_consts[104], mod_consts[108], NULL, 0, 0, 0);
}

// The module function declarations.
static PyObject *MAKE_FUNCTION___main__$$$function__1_main();


// The module function definitions.
static PyObject *impl___main__$$$function__1_main(struct Nuitka_FunctionObject const *self, PyObject **python_pars) {
    // Preserve error status for checks
#ifndef __NUITKA_NO_ASSERT__
    NUITKA_MAY_BE_UNUSED bool had_error = ERROR_OCCURRED();
#endif

    // Local variable declarations.
    PyObject *var_parser = NULL;
    PyObject *var_requiredNamed = NULL;
    PyObject *var_args = NULL;
    PyObject *var_plasmid = NULL;
    PyObject *var_rank = NULL;
    PyObject *var_tree = NULL;
    PyObject *var_Prediction_object = NULL;
    PyObject *var_Phylogenetic_object = NULL;
    PyObject *var_input_file = NULL;
    PyObject *var_files_lines = NULL;
    PyObject *var_unique_sequences = NULL;
    PyObject *var_input_sequence_list = NULL;
    PyObject *var_label_list = NULL;
    PyObject *var_prediction_output = NULL;
    PyObject *var_significance_output = NULL;
    PyObject *var_taxa_output = NULL;
    PyObject *var_converted_sequences_phyolgenetic_tree = NULL;
    PyObject *var_j = NULL;
    PyObject *var_bar = NULL;
    PyObject *var_user = NULL;
    PyObject *var_unique_kmer_frequency = NULL;
    PyObject *var_distance_matrix = NULL;
    PyObject *outline_0_var_i = NULL;
    PyObject *tmp_for_loop_1__for_iterator = NULL;
    PyObject *tmp_for_loop_1__iter_value = NULL;
    PyObject *tmp_for_loop_2__for_iterator = NULL;
    PyObject *tmp_for_loop_2__iter_value = NULL;
    PyObject *tmp_listcomp_1__$0 = NULL;
    PyObject *tmp_listcomp_1__contraction = NULL;
    PyObject *tmp_listcomp_1__iter_value_0 = NULL;
    PyObject *tmp_with_1__enter = NULL;
    PyObject *tmp_with_1__exit = NULL;
    nuitka_bool tmp_with_1__indicator = NUITKA_BOOL_UNASSIGNED;
    PyObject *tmp_with_1__source = NULL;
    struct Nuitka_FrameObject *frame_7fdb2c580e4cfeb4386220f78b54ec88;
    NUITKA_MAY_BE_UNUSED char const *type_description_1 = NULL;
    PyObject *exception_type = NULL;
    PyObject *exception_value = NULL;
    PyTracebackObject *exception_tb = NULL;
    NUITKA_MAY_BE_UNUSED int exception_lineno = 0;
    NUITKA_MAY_BE_UNUSED nuitka_void tmp_unused;
    PyObject *exception_keeper_type_1;
    PyObject *exception_keeper_value_1;
    PyTracebackObject *exception_keeper_tb_1;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_1;
    PyObject *exception_keeper_type_2;
    PyObject *exception_keeper_value_2;
    PyTracebackObject *exception_keeper_tb_2;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_2;
    PyObject *exception_preserved_type_1;
    PyObject *exception_preserved_value_1;
    PyTracebackObject *exception_preserved_tb_1;
    int tmp_res;
    bool tmp_result;
    PyObject *exception_keeper_type_3;
    PyObject *exception_keeper_value_3;
    PyTracebackObject *exception_keeper_tb_3;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_3;
    PyObject *exception_keeper_type_4;
    PyObject *exception_keeper_value_4;
    PyTracebackObject *exception_keeper_tb_4;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_4;
    PyObject *exception_keeper_type_5;
    PyObject *exception_keeper_value_5;
    PyTracebackObject *exception_keeper_tb_5;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_5;
    struct Nuitka_FrameObject *frame_ccd1fb807920dfdbbfd5c88e344983c6_2;
    NUITKA_MAY_BE_UNUSED char const *type_description_2 = NULL;
    PyObject *exception_keeper_type_6;
    PyObject *exception_keeper_value_6;
    PyTracebackObject *exception_keeper_tb_6;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_6;
    static struct Nuitka_FrameObject *cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2 = NULL;
    PyObject *exception_keeper_type_7;
    PyObject *exception_keeper_value_7;
    PyTracebackObject *exception_keeper_tb_7;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_7;
    PyObject *exception_keeper_type_8;
    PyObject *exception_keeper_value_8;
    PyTracebackObject *exception_keeper_tb_8;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_8;
    static struct Nuitka_FrameObject *cache_frame_7fdb2c580e4cfeb4386220f78b54ec88 = NULL;
    PyObject *tmp_return_value = NULL;
    PyObject *exception_keeper_type_9;
    PyObject *exception_keeper_value_9;
    PyTracebackObject *exception_keeper_tb_9;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_9;

    // Actual function body.
    // Tried code:
    if (isFrameUnusable(cache_frame_7fdb2c580e4cfeb4386220f78b54ec88)) {
        Py_XDECREF(cache_frame_7fdb2c580e4cfeb4386220f78b54ec88);

#if _DEBUG_REFCOUNTS
        if (cache_frame_7fdb2c580e4cfeb4386220f78b54ec88 == NULL) {
            count_active_frame_cache_instances += 1;
        } else {
            count_released_frame_cache_instances += 1;
        }
        count_allocated_frame_cache_instances += 1;
#endif
        cache_frame_7fdb2c580e4cfeb4386220f78b54ec88 = MAKE_FUNCTION_FRAME(codeobj_7fdb2c580e4cfeb4386220f78b54ec88, module___main__, sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *)+sizeof(void *));
#if _DEBUG_REFCOUNTS
    } else {
        count_hit_frame_cache_instances += 1;
#endif
    }
    assert(cache_frame_7fdb2c580e4cfeb4386220f78b54ec88->m_type_description == NULL);
    frame_7fdb2c580e4cfeb4386220f78b54ec88 = cache_frame_7fdb2c580e4cfeb4386220f78b54ec88;

    // Push the new frame as the currently active one.
    pushFrameStack(frame_7fdb2c580e4cfeb4386220f78b54ec88);

    // Mark the frame object as in use, ref count 1 will be up for reuse.
    assert(Py_REFCNT(frame_7fdb2c580e4cfeb4386220f78b54ec88) == 2); // Frame stack

    // Framed code:
    {
        PyObject *tmp_assign_source_1;
        PyObject *tmp_called_value_1;
        PyObject *tmp_kw_call_value_0_1;
        tmp_called_value_1 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[0]);

        if (unlikely(tmp_called_value_1 == NULL)) {
            tmp_called_value_1 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[0]);
        }

        if (tmp_called_value_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 25;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_kw_call_value_0_1 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[1]);

        if (unlikely(tmp_kw_call_value_0_1 == NULL)) {
            tmp_kw_call_value_0_1 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[1]);
        }

        if (tmp_kw_call_value_0_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 25;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 25;
        {
            PyObject *kw_values[1] = {tmp_kw_call_value_0_1};

            tmp_assign_source_1 = CALL_FUNCTION_WITH_NO_ARGS_KWSPLIT(tmp_called_value_1, kw_values, mod_consts[2]);
        }

        if (tmp_assign_source_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 25;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_parser == NULL);
        var_parser = tmp_assign_source_1;
    }
    {
        PyObject *tmp_assign_source_2;
        PyObject *tmp_called_instance_1;
        CHECK_OBJECT(var_parser);
        tmp_called_instance_1 = var_parser;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 26;
        tmp_assign_source_2 = CALL_METHOD_WITH_SINGLE_ARG(
            tmp_called_instance_1,
            mod_consts[3],
            PyTuple_GET_ITEM(mod_consts[4], 0)
        );

        if (tmp_assign_source_2 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 26;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_requiredNamed == NULL);
        var_requiredNamed = tmp_assign_source_2;
    }
    {
        PyObject *tmp_called_value_2;
        PyObject *tmp_expression_value_1;
        PyObject *tmp_call_result_1;
        CHECK_OBJECT(var_requiredNamed);
        tmp_expression_value_1 = var_requiredNamed;
        tmp_called_value_2 = LOOKUP_ATTRIBUTE(tmp_expression_value_1, mod_consts[5]);
        if (tmp_called_value_2 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 27;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 27;
        tmp_call_result_1 = CALL_FUNCTION_WITH_ARGS2_VECTORCALL(tmp_called_value_2, &PyTuple_GET_ITEM(mod_consts[6], 0), mod_consts[7]);
        Py_DECREF(tmp_called_value_2);
        if (tmp_call_result_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 27;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_1);
    }
    {
        PyObject *tmp_called_value_3;
        PyObject *tmp_expression_value_2;
        PyObject *tmp_call_result_2;
        CHECK_OBJECT(var_requiredNamed);
        tmp_expression_value_2 = var_requiredNamed;
        tmp_called_value_3 = LOOKUP_ATTRIBUTE(tmp_expression_value_2, mod_consts[5]);
        if (tmp_called_value_3 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 28;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 28;
        tmp_call_result_2 = CALL_FUNCTION_WITH_ARGS2_VECTORCALL(tmp_called_value_3, &PyTuple_GET_ITEM(mod_consts[8], 0), mod_consts[9]);
        Py_DECREF(tmp_called_value_3);
        if (tmp_call_result_2 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 28;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_2);
    }
    {
        PyObject *tmp_called_value_4;
        PyObject *tmp_expression_value_3;
        PyObject *tmp_call_result_3;
        PyObject *tmp_call_args_values_1;
        CHECK_OBJECT(var_parser);
        tmp_expression_value_3 = var_parser;
        tmp_called_value_4 = LOOKUP_ATTRIBUTE(tmp_expression_value_3, mod_consts[5]);
        if (tmp_called_value_4 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 30;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_call_args_values_1 = DEEP_COPY_TUPLE(mod_consts[10]);
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 30;
        tmp_call_result_3 = CALL_FUNCTION_WITH_ARGS2_VECTORCALL(tmp_called_value_4, &PyTuple_GET_ITEM(tmp_call_args_values_1, 0), mod_consts[11]);
        Py_DECREF(tmp_called_value_4);
        Py_DECREF(tmp_call_args_values_1);
        if (tmp_call_result_3 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 30;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_3);
    }
    {
        PyObject *tmp_called_value_5;
        PyObject *tmp_expression_value_4;
        PyObject *tmp_call_result_4;
        PyObject *tmp_call_args_values_2;
        CHECK_OBJECT(var_parser);
        tmp_expression_value_4 = var_parser;
        tmp_called_value_5 = LOOKUP_ATTRIBUTE(tmp_expression_value_4, mod_consts[5]);
        if (tmp_called_value_5 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 31;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_call_args_values_2 = DEEP_COPY_TUPLE(mod_consts[12]);
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 31;
        tmp_call_result_4 = CALL_FUNCTION_WITH_ARGS2_VECTORCALL(tmp_called_value_5, &PyTuple_GET_ITEM(tmp_call_args_values_2, 0), mod_consts[11]);
        Py_DECREF(tmp_called_value_5);
        Py_DECREF(tmp_call_args_values_2);
        if (tmp_call_result_4 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 31;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_4);
    }
    {
        PyObject *tmp_called_value_6;
        PyObject *tmp_expression_value_5;
        PyObject *tmp_call_result_5;
        CHECK_OBJECT(var_parser);
        tmp_expression_value_5 = var_parser;
        tmp_called_value_6 = LOOKUP_ATTRIBUTE(tmp_expression_value_5, mod_consts[5]);
        if (tmp_called_value_6 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 32;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 32;
        tmp_call_result_5 = CALL_FUNCTION_WITH_ARGS2_VECTORCALL(tmp_called_value_6, &PyTuple_GET_ITEM(mod_consts[13], 0), mod_consts[14]);
        Py_DECREF(tmp_called_value_6);
        if (tmp_call_result_5 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 32;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_5);
    }
    {
        PyObject *tmp_called_value_7;
        PyObject *tmp_expression_value_6;
        PyObject *tmp_call_result_6;
        CHECK_OBJECT(var_parser);
        tmp_expression_value_6 = var_parser;
        tmp_called_value_7 = LOOKUP_ATTRIBUTE(tmp_expression_value_6, mod_consts[15]);
        if (tmp_called_value_7 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 33;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 33;
        tmp_call_result_6 = CALL_FUNCTION_WITH_NO_ARGS_KWSPLIT(tmp_called_value_7, &PyTuple_GET_ITEM(mod_consts[16], 0), mod_consts[17]);
        Py_DECREF(tmp_called_value_7);
        if (tmp_call_result_6 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 33;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_6);
    }
    {
        PyObject *tmp_assign_source_3;
        PyObject *tmp_vars_arg_1;
        PyObject *tmp_called_instance_2;
        CHECK_OBJECT(var_parser);
        tmp_called_instance_2 = var_parser;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 34;
        tmp_vars_arg_1 = CALL_METHOD_NO_ARGS(tmp_called_instance_2, mod_consts[18]);
        if (tmp_vars_arg_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 34;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_assign_source_3 = LOOKUP_VARS(tmp_vars_arg_1);
        Py_DECREF(tmp_vars_arg_1);
        if (tmp_assign_source_3 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 34;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_args == NULL);
        var_args = tmp_assign_source_3;
    }
    {
        PyObject *tmp_assign_source_4;
        PyObject *tmp_expression_value_7;
        PyObject *tmp_subscript_value_1;
        CHECK_OBJECT(var_args);
        tmp_expression_value_7 = var_args;
        tmp_subscript_value_1 = mod_consts[19];
        tmp_assign_source_4 = LOOKUP_SUBSCRIPT(tmp_expression_value_7, tmp_subscript_value_1);
        if (tmp_assign_source_4 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 36;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_plasmid == NULL);
        var_plasmid = tmp_assign_source_4;
    }
    {
        PyObject *tmp_assign_source_5;
        PyObject *tmp_expression_value_8;
        PyObject *tmp_subscript_value_2;
        CHECK_OBJECT(var_args);
        tmp_expression_value_8 = var_args;
        tmp_subscript_value_2 = mod_consts[20];
        tmp_assign_source_5 = LOOKUP_SUBSCRIPT(tmp_expression_value_8, tmp_subscript_value_2);
        if (tmp_assign_source_5 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 37;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_rank == NULL);
        var_rank = tmp_assign_source_5;
    }
    {
        PyObject *tmp_assign_source_6;
        PyObject *tmp_expression_value_9;
        PyObject *tmp_subscript_value_3;
        CHECK_OBJECT(var_args);
        tmp_expression_value_9 = var_args;
        tmp_subscript_value_3 = mod_consts[21];
        tmp_assign_source_6 = LOOKUP_SUBSCRIPT(tmp_expression_value_9, tmp_subscript_value_3);
        if (tmp_assign_source_6 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 38;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_tree == NULL);
        var_tree = tmp_assign_source_6;
    }
    {
        PyObject *tmp_assign_source_7;
        PyObject *tmp_called_value_8;
        PyObject *tmp_args_element_value_1;
        tmp_called_value_8 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[22]);

        if (unlikely(tmp_called_value_8 == NULL)) {
            tmp_called_value_8 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[22]);
        }

        if (tmp_called_value_8 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 41;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_args_element_value_1 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[23]);

        if (unlikely(tmp_args_element_value_1 == NULL)) {
            tmp_args_element_value_1 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[23]);
        }

        if (tmp_args_element_value_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 41;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 41;
        tmp_assign_source_7 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_8, tmp_args_element_value_1);
        if (tmp_assign_source_7 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 41;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_Prediction_object == NULL);
        var_Prediction_object = tmp_assign_source_7;
    }
    {
        PyObject *tmp_assign_source_8;
        tmp_assign_source_8 = PyList_New(0);
        assert(var_Phylogenetic_object == NULL);
        var_Phylogenetic_object = tmp_assign_source_8;
    }
    {
        nuitka_bool tmp_condition_result_1;
        int tmp_truth_name_1;
        CHECK_OBJECT(var_tree);
        tmp_truth_name_1 = CHECK_IF_TRUE(var_tree);
        if (tmp_truth_name_1 == -1) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 43;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_condition_result_1 = tmp_truth_name_1 == 0 ? NUITKA_BOOL_FALSE : NUITKA_BOOL_TRUE;
        if (tmp_condition_result_1 == NUITKA_BOOL_TRUE) {
            goto branch_yes_1;
        } else {
            goto branch_no_1;
        }
    }
    branch_yes_1:;
    {
        PyObject *tmp_assign_source_9;
        PyObject *tmp_called_value_9;
        PyObject *tmp_args_element_value_2;
        tmp_called_value_9 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[24]);

        if (unlikely(tmp_called_value_9 == NULL)) {
            tmp_called_value_9 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[24]);
        }

        if (tmp_called_value_9 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 44;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_args_element_value_2 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[23]);

        if (unlikely(tmp_args_element_value_2 == NULL)) {
            tmp_args_element_value_2 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[23]);
        }

        if (tmp_args_element_value_2 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 44;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 44;
        tmp_assign_source_9 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_9, tmp_args_element_value_2);
        if (tmp_assign_source_9 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 44;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        {
            PyObject *old = var_Phylogenetic_object;
            assert(old != NULL);
            var_Phylogenetic_object = tmp_assign_source_9;
            Py_DECREF(old);
        }

    }
    branch_no_1:;
    {
        PyObject *tmp_called_value_10;
        PyObject *tmp_call_result_7;
        tmp_called_value_10 = LOOKUP_BUILTIN(mod_consts[25]);
        assert(tmp_called_value_10 != NULL);
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 46;
        tmp_call_result_7 = CALL_FUNCTION_WITH_POSARGS1(tmp_called_value_10, mod_consts[26]);

        if (tmp_call_result_7 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 46;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_7);
    }
    {
        PyObject *tmp_assign_source_10;
        PyObject *tmp_open_filename_1;
        PyObject *tmp_called_value_11;
        PyObject *tmp_expression_value_10;
        PyObject *tmp_expression_value_11;
        PyObject *tmp_args_element_value_3;
        PyObject *tmp_args_element_value_4;
        PyObject *tmp_expression_value_12;
        PyObject *tmp_subscript_value_4;
        PyObject *tmp_open_mode_1;
        tmp_expression_value_11 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[27]);

        if (unlikely(tmp_expression_value_11 == NULL)) {
            tmp_expression_value_11 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[27]);
        }

        assert(!(tmp_expression_value_11 == NULL));
        tmp_expression_value_10 = LOOKUP_ATTRIBUTE(tmp_expression_value_11, mod_consts[28]);
        if (tmp_expression_value_10 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 47;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_called_value_11 = LOOKUP_ATTRIBUTE(tmp_expression_value_10, mod_consts[29]);
        Py_DECREF(tmp_expression_value_10);
        if (tmp_called_value_11 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 47;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_args_element_value_3 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[23]);

        if (unlikely(tmp_args_element_value_3 == NULL)) {
            tmp_args_element_value_3 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[23]);
        }

        if (tmp_args_element_value_3 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_11);

            exception_lineno = 47;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        CHECK_OBJECT(var_args);
        tmp_expression_value_12 = var_args;
        tmp_subscript_value_4 = mod_consts[30];
        tmp_args_element_value_4 = LOOKUP_SUBSCRIPT(tmp_expression_value_12, tmp_subscript_value_4);
        if (tmp_args_element_value_4 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_11);

            exception_lineno = 47;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 47;
        {
            PyObject *call_args[] = {tmp_args_element_value_3, tmp_args_element_value_4};
            tmp_open_filename_1 = CALL_FUNCTION_WITH_ARGS2(tmp_called_value_11, call_args);
        }

        Py_DECREF(tmp_called_value_11);
        Py_DECREF(tmp_args_element_value_4);
        if (tmp_open_filename_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 47;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_open_mode_1 = mod_consts[31];
        tmp_assign_source_10 = BUILTIN_OPEN(tmp_open_filename_1, tmp_open_mode_1, NULL, NULL, NULL, NULL, NULL, NULL);
        Py_DECREF(tmp_open_filename_1);
        if (tmp_assign_source_10 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 47;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_input_file == NULL);
        var_input_file = tmp_assign_source_10;
    }
    {
        PyObject *tmp_assign_source_11;
        PyObject *tmp_called_instance_3;
        CHECK_OBJECT(var_input_file);
        tmp_called_instance_3 = var_input_file;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 48;
        tmp_assign_source_11 = CALL_METHOD_NO_ARGS(tmp_called_instance_3, mod_consts[32]);
        if (tmp_assign_source_11 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 48;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_files_lines == NULL);
        var_files_lines = tmp_assign_source_11;
    }
    {
        PyObject *tmp_assign_source_12;
        PyObject *tmp_called_value_12;
        PyObject *tmp_expression_value_13;
        PyObject *tmp_args_element_value_5;
        PyObject *tmp_called_value_13;
        PyObject *tmp_expression_value_14;
        PyObject *tmp_expression_value_15;
        PyObject *tmp_args_element_value_6;
        PyObject *tmp_args_element_value_7;
        PyObject *tmp_args_element_value_8;
        PyObject *tmp_left_value_1;
        PyObject *tmp_left_value_2;
        PyObject *tmp_right_value_1;
        PyObject *tmp_right_value_2;
        tmp_expression_value_13 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[33]);

        if (unlikely(tmp_expression_value_13 == NULL)) {
            tmp_expression_value_13 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[33]);
        }

        if (tmp_expression_value_13 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_called_value_12 = LOOKUP_ATTRIBUTE(tmp_expression_value_13, mod_consts[34]);
        if (tmp_called_value_12 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_expression_value_15 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[27]);

        if (unlikely(tmp_expression_value_15 == NULL)) {
            tmp_expression_value_15 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[27]);
        }

        if (tmp_expression_value_15 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_12);

            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_expression_value_14 = LOOKUP_ATTRIBUTE(tmp_expression_value_15, mod_consts[28]);
        if (tmp_expression_value_14 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_12);

            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_called_value_13 = LOOKUP_ATTRIBUTE(tmp_expression_value_14, mod_consts[29]);
        Py_DECREF(tmp_expression_value_14);
        if (tmp_called_value_13 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_12);

            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_args_element_value_6 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[23]);

        if (unlikely(tmp_args_element_value_6 == NULL)) {
            tmp_args_element_value_6 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[23]);
        }

        if (tmp_args_element_value_6 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_12);
            Py_DECREF(tmp_called_value_13);

            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_args_element_value_7 = mod_consts[35];
        tmp_left_value_2 = mod_consts[36];
        CHECK_OBJECT(var_plasmid);
        tmp_right_value_1 = var_plasmid;
        tmp_left_value_1 = BINARY_OPERATION_ADD_OBJECT_UNICODE_OBJECT(tmp_left_value_2, tmp_right_value_1);
        if (tmp_left_value_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_12);
            Py_DECREF(tmp_called_value_13);

            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_right_value_2 = mod_consts[37];
        tmp_args_element_value_8 = BINARY_OPERATION_ADD_OBJECT_OBJECT_UNICODE(tmp_left_value_1, tmp_right_value_2);
        Py_DECREF(tmp_left_value_1);
        if (tmp_args_element_value_8 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_12);
            Py_DECREF(tmp_called_value_13);

            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 51;
        {
            PyObject *call_args[] = {tmp_args_element_value_6, tmp_args_element_value_7, tmp_args_element_value_8};
            tmp_args_element_value_5 = CALL_FUNCTION_WITH_ARGS3(tmp_called_value_13, call_args);
        }

        Py_DECREF(tmp_called_value_13);
        Py_DECREF(tmp_args_element_value_8);
        if (tmp_args_element_value_5 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_12);

            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 51;
        tmp_assign_source_12 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_12, tmp_args_element_value_5);
        Py_DECREF(tmp_called_value_12);
        Py_DECREF(tmp_args_element_value_5);
        if (tmp_assign_source_12 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 51;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_unique_sequences == NULL);
        var_unique_sequences = tmp_assign_source_12;
    }
    {
        PyObject *tmp_assign_source_13;
        PyObject *tmp_expression_value_16;
        PyObject *tmp_expression_value_17;
        PyObject *tmp_expression_value_18;
        PyObject *tmp_subscript_value_5;
        CHECK_OBJECT(var_unique_sequences);
        tmp_expression_value_18 = var_unique_sequences;
        tmp_expression_value_17 = LOOKUP_ATTRIBUTE(tmp_expression_value_18, mod_consts[38]);
        if (tmp_expression_value_17 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 52;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_subscript_value_5 = mod_consts[39];
        tmp_expression_value_16 = LOOKUP_SUBSCRIPT(tmp_expression_value_17, tmp_subscript_value_5);
        Py_DECREF(tmp_expression_value_17);
        if (tmp_expression_value_16 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 52;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_assign_source_13 = LOOKUP_ATTRIBUTE(tmp_expression_value_16, mod_consts[40]);
        Py_DECREF(tmp_expression_value_16);
        if (tmp_assign_source_13 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 52;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        {
            PyObject *old = var_unique_sequences;
            assert(old != NULL);
            var_unique_sequences = tmp_assign_source_13;
            Py_DECREF(old);
        }

    }
    {
        PyObject *tmp_assign_source_14;
        PyObject *tmp_called_value_14;
        PyObject *tmp_expression_value_19;
        PyObject *tmp_args_value_1;
        PyObject *tmp_tuple_element_1;
        PyObject *tmp_kwargs_value_1;
        CHECK_OBJECT(var_Prediction_object);
        tmp_expression_value_19 = var_Prediction_object;
        tmp_called_value_14 = LOOKUP_ATTRIBUTE(tmp_expression_value_19, mod_consts[41]);
        if (tmp_called_value_14 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 54;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        CHECK_OBJECT(var_files_lines);
        tmp_tuple_element_1 = var_files_lines;
        tmp_args_value_1 = PyTuple_New(1);
        PyTuple_SET_ITEM0(tmp_args_value_1, 0, tmp_tuple_element_1);
        tmp_kwargs_value_1 = PyDict_Copy(mod_consts[42]);
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 54;
        tmp_assign_source_14 = CALL_FUNCTION(tmp_called_value_14, tmp_args_value_1, tmp_kwargs_value_1);
        Py_DECREF(tmp_called_value_14);
        Py_DECREF(tmp_args_value_1);
        Py_DECREF(tmp_kwargs_value_1);
        if (tmp_assign_source_14 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 54;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_input_sequence_list == NULL);
        var_input_sequence_list = tmp_assign_source_14;
    }
    {
        PyObject *tmp_assign_source_15;
        PyObject *tmp_called_value_15;
        PyObject *tmp_expression_value_20;
        PyObject *tmp_args_value_2;
        PyObject *tmp_tuple_element_2;
        PyObject *tmp_kwargs_value_2;
        CHECK_OBJECT(var_Prediction_object);
        tmp_expression_value_20 = var_Prediction_object;
        tmp_called_value_15 = LOOKUP_ATTRIBUTE(tmp_expression_value_20, mod_consts[41]);
        if (tmp_called_value_15 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 55;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        CHECK_OBJECT(var_files_lines);
        tmp_tuple_element_2 = var_files_lines;
        tmp_args_value_2 = PyTuple_New(1);
        PyTuple_SET_ITEM0(tmp_args_value_2, 0, tmp_tuple_element_2);
        tmp_kwargs_value_2 = PyDict_Copy(mod_consts[43]);
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 55;
        tmp_assign_source_15 = CALL_FUNCTION(tmp_called_value_15, tmp_args_value_2, tmp_kwargs_value_2);
        Py_DECREF(tmp_called_value_15);
        Py_DECREF(tmp_args_value_2);
        Py_DECREF(tmp_kwargs_value_2);
        if (tmp_assign_source_15 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 55;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_label_list == NULL);
        var_label_list = tmp_assign_source_15;
    }
    {
        PyObject *tmp_assign_source_16;
        tmp_assign_source_16 = PyList_New(0);
        assert(var_prediction_output == NULL);
        var_prediction_output = tmp_assign_source_16;
    }
    {
        PyObject *tmp_assign_source_17;
        tmp_assign_source_17 = PyList_New(0);
        assert(var_significance_output == NULL);
        var_significance_output = tmp_assign_source_17;
    }
    {
        PyObject *tmp_assign_source_18;
        tmp_assign_source_18 = PyList_New(0);
        assert(var_taxa_output == NULL);
        var_taxa_output = tmp_assign_source_18;
    }
    {
        PyObject *tmp_assign_source_19;
        tmp_assign_source_19 = PyList_New(0);
        assert(var_converted_sequences_phyolgenetic_tree == NULL);
        var_converted_sequences_phyolgenetic_tree = tmp_assign_source_19;
    }
    {
        PyObject *tmp_called_value_16;
        PyObject *tmp_call_result_8;
        PyObject *tmp_args_element_value_9;
        PyObject *tmp_left_value_3;
        PyObject *tmp_right_value_3;
        tmp_called_value_16 = LOOKUP_BUILTIN(mod_consts[25]);
        assert(tmp_called_value_16 != NULL);
        tmp_left_value_3 = mod_consts[44];
        CHECK_OBJECT(var_plasmid);
        tmp_right_value_3 = var_plasmid;
        tmp_args_element_value_9 = BINARY_OPERATION_ADD_OBJECT_UNICODE_OBJECT(tmp_left_value_3, tmp_right_value_3);
        if (tmp_args_element_value_9 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 62;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 62;
        tmp_call_result_8 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_16, tmp_args_element_value_9);
        Py_DECREF(tmp_args_element_value_9);
        if (tmp_call_result_8 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 62;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_8);
    }
    {
        PyObject *tmp_assign_source_20;
        PyObject *tmp_iter_arg_1;
        PyObject *tmp_xrange_low_1;
        PyObject *tmp_len_arg_1;
        CHECK_OBJECT(var_input_sequence_list);
        tmp_len_arg_1 = var_input_sequence_list;
        tmp_xrange_low_1 = BUILTIN_LEN(tmp_len_arg_1);
        if (tmp_xrange_low_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 63;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_iter_arg_1 = BUILTIN_XRANGE1(tmp_xrange_low_1);
        Py_DECREF(tmp_xrange_low_1);
        if (tmp_iter_arg_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 63;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_assign_source_20 = MAKE_ITERATOR(tmp_iter_arg_1);
        Py_DECREF(tmp_iter_arg_1);
        if (tmp_assign_source_20 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 63;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(tmp_for_loop_1__for_iterator == NULL);
        tmp_for_loop_1__for_iterator = tmp_assign_source_20;
    }
    // Tried code:
    loop_start_1:;
    {
        PyObject *tmp_next_source_1;
        PyObject *tmp_assign_source_21;
        CHECK_OBJECT(tmp_for_loop_1__for_iterator);
        tmp_next_source_1 = tmp_for_loop_1__for_iterator;
        tmp_assign_source_21 = ITERATOR_NEXT(tmp_next_source_1);
        if (tmp_assign_source_21 == NULL) {
            if (CHECK_AND_CLEAR_STOP_ITERATION_OCCURRED()) {

                goto loop_end_1;
            } else {

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
                type_description_1 = "oooooooooooooooooooooo";
                exception_lineno = 63;
                goto try_except_handler_2;
            }
        }

        {
            PyObject *old = tmp_for_loop_1__iter_value;
            tmp_for_loop_1__iter_value = tmp_assign_source_21;
            Py_XDECREF(old);
        }

    }
    {
        PyObject *tmp_assign_source_22;
        CHECK_OBJECT(tmp_for_loop_1__iter_value);
        tmp_assign_source_22 = tmp_for_loop_1__iter_value;
        {
            PyObject *old = var_j;
            var_j = tmp_assign_source_22;
            Py_INCREF(var_j);
            Py_XDECREF(old);
        }

    }
    {
        PyObject *tmp_called_value_17;
        PyObject *tmp_call_result_9;
        PyObject *tmp_args_element_value_10;
        PyObject *tmp_left_value_4;
        PyObject *tmp_right_value_4;
        PyObject *tmp_expression_value_21;
        PyObject *tmp_subscript_value_6;
        tmp_called_value_17 = LOOKUP_BUILTIN(mod_consts[25]);
        assert(tmp_called_value_17 != NULL);
        tmp_left_value_4 = mod_consts[45];
        CHECK_OBJECT(var_label_list);
        tmp_expression_value_21 = var_label_list;
        CHECK_OBJECT(var_j);
        tmp_subscript_value_6 = var_j;
        tmp_right_value_4 = LOOKUP_SUBSCRIPT(tmp_expression_value_21, tmp_subscript_value_6);
        if (tmp_right_value_4 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 64;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        tmp_args_element_value_10 = BINARY_OPERATION_ADD_OBJECT_UNICODE_OBJECT(tmp_left_value_4, tmp_right_value_4);
        Py_DECREF(tmp_right_value_4);
        if (tmp_args_element_value_10 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 64;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 64;
        tmp_call_result_9 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_17, tmp_args_element_value_10);
        Py_DECREF(tmp_args_element_value_10);
        if (tmp_call_result_9 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 64;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        Py_DECREF(tmp_call_result_9);
    }
    // Tried code:
    {
        PyObject *tmp_assign_source_23;
        PyObject *tmp_called_value_18;
        PyObject *tmp_expression_value_22;
        tmp_expression_value_22 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[46]);

        if (unlikely(tmp_expression_value_22 == NULL)) {
            tmp_expression_value_22 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[46]);
        }

        if (tmp_expression_value_22 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 65;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_3;
        }
        tmp_called_value_18 = LOOKUP_ATTRIBUTE(tmp_expression_value_22, mod_consts[47]);
        if (tmp_called_value_18 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 65;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_3;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 65;
        tmp_assign_source_23 = CALL_FUNCTION_WITH_ARGS1_VECTORCALL(tmp_called_value_18, &PyTuple_GET_ITEM(mod_consts[48], 0), mod_consts[49]);
        Py_DECREF(tmp_called_value_18);
        if (tmp_assign_source_23 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 65;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_3;
        }
        {
            PyObject *old = tmp_with_1__source;
            tmp_with_1__source = tmp_assign_source_23;
            Py_XDECREF(old);
        }

    }
    {
        PyObject *tmp_assign_source_24;
        PyObject *tmp_called_value_19;
        PyObject *tmp_expression_value_23;
        CHECK_OBJECT(tmp_with_1__source);
        tmp_expression_value_23 = tmp_with_1__source;
        tmp_called_value_19 = LOOKUP_SPECIAL(tmp_expression_value_23, mod_consts[50]);
        if (tmp_called_value_19 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 65;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_3;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 65;
        tmp_assign_source_24 = CALL_FUNCTION_NO_ARGS(tmp_called_value_19);
        Py_DECREF(tmp_called_value_19);
        if (tmp_assign_source_24 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 65;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_3;
        }
        {
            PyObject *old = tmp_with_1__enter;
            tmp_with_1__enter = tmp_assign_source_24;
            Py_XDECREF(old);
        }

    }
    {
        PyObject *tmp_assign_source_25;
        PyObject *tmp_expression_value_24;
        CHECK_OBJECT(tmp_with_1__source);
        tmp_expression_value_24 = tmp_with_1__source;
        tmp_assign_source_25 = LOOKUP_SPECIAL(tmp_expression_value_24, mod_consts[51]);
        if (tmp_assign_source_25 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 65;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_3;
        }
        {
            PyObject *old = tmp_with_1__exit;
            tmp_with_1__exit = tmp_assign_source_25;
            Py_XDECREF(old);
        }

    }
    {
        nuitka_bool tmp_assign_source_26;
        tmp_assign_source_26 = NUITKA_BOOL_TRUE;
        tmp_with_1__indicator = tmp_assign_source_26;
    }
    {
        PyObject *tmp_assign_source_27;
        CHECK_OBJECT(tmp_with_1__enter);
        tmp_assign_source_27 = tmp_with_1__enter;
        {
            PyObject *old = var_bar;
            var_bar = tmp_assign_source_27;
            Py_INCREF(var_bar);
            Py_XDECREF(old);
        }

    }
    // Tried code:
    // Tried code:
    {
        PyObject *tmp_assign_source_28;
        PyObject *tmp_iter_arg_2;
        CHECK_OBJECT(var_bar);
        tmp_iter_arg_2 = var_bar;
        tmp_assign_source_28 = MAKE_ITERATOR(tmp_iter_arg_2);
        if (tmp_assign_source_28 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 66;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_5;
        }
        {
            PyObject *old = tmp_for_loop_2__for_iterator;
            tmp_for_loop_2__for_iterator = tmp_assign_source_28;
            Py_XDECREF(old);
        }

    }
    // Tried code:
    loop_start_2:;
    {
        PyObject *tmp_next_source_2;
        PyObject *tmp_assign_source_29;
        CHECK_OBJECT(tmp_for_loop_2__for_iterator);
        tmp_next_source_2 = tmp_for_loop_2__for_iterator;
        tmp_assign_source_29 = ITERATOR_NEXT(tmp_next_source_2);
        if (tmp_assign_source_29 == NULL) {
            if (CHECK_AND_CLEAR_STOP_ITERATION_OCCURRED()) {

                goto loop_end_2;
            } else {

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
                type_description_1 = "oooooooooooooooooooooo";
                exception_lineno = 66;
                goto try_except_handler_6;
            }
        }

        {
            PyObject *old = tmp_for_loop_2__iter_value;
            tmp_for_loop_2__iter_value = tmp_assign_source_29;
            Py_XDECREF(old);
        }

    }
    {
        PyObject *tmp_assign_source_30;
        CHECK_OBJECT(tmp_for_loop_2__iter_value);
        tmp_assign_source_30 = tmp_for_loop_2__iter_value;
        {
            PyObject *old = var_user;
            var_user = tmp_assign_source_30;
            Py_INCREF(var_user);
            Py_XDECREF(old);
        }

    }
    {
        PyObject *tmp_called_instance_4;
        PyObject *tmp_call_result_10;
        tmp_called_instance_4 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[52]);

        if (unlikely(tmp_called_instance_4 == NULL)) {
            tmp_called_instance_4 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[52]);
        }

        if (tmp_called_instance_4 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 67;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_6;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 67;
        tmp_call_result_10 = CALL_METHOD_WITH_SINGLE_ARG(
            tmp_called_instance_4,
            mod_consts[53],
            PyTuple_GET_ITEM(mod_consts[54], 0)
        );

        if (tmp_call_result_10 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 67;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_6;
        }
        Py_DECREF(tmp_call_result_10);
    }
    if (CONSIDER_THREADING() == false) {
        assert(ERROR_OCCURRED());

        FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


        exception_lineno = 66;
        type_description_1 = "oooooooooooooooooooooo";
        goto try_except_handler_6;
    }
    goto loop_start_2;
    loop_end_2:;
    goto try_end_1;
    // Exception handler code:
    try_except_handler_6:;
    exception_keeper_type_1 = exception_type;
    exception_keeper_value_1 = exception_value;
    exception_keeper_tb_1 = exception_tb;
    exception_keeper_lineno_1 = exception_lineno;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    exception_lineno = 0;

    Py_XDECREF(tmp_for_loop_2__iter_value);
    tmp_for_loop_2__iter_value = NULL;
    CHECK_OBJECT(tmp_for_loop_2__for_iterator);
    Py_DECREF(tmp_for_loop_2__for_iterator);
    tmp_for_loop_2__for_iterator = NULL;
    // Re-raise.
    exception_type = exception_keeper_type_1;
    exception_value = exception_keeper_value_1;
    exception_tb = exception_keeper_tb_1;
    exception_lineno = exception_keeper_lineno_1;

    goto try_except_handler_5;
    // End of try:
    try_end_1:;
    Py_XDECREF(tmp_for_loop_2__iter_value);
    tmp_for_loop_2__iter_value = NULL;
    CHECK_OBJECT(tmp_for_loop_2__for_iterator);
    Py_DECREF(tmp_for_loop_2__for_iterator);
    tmp_for_loop_2__for_iterator = NULL;
    goto try_end_2;
    // Exception handler code:
    try_except_handler_5:;
    exception_keeper_type_2 = exception_type;
    exception_keeper_value_2 = exception_value;
    exception_keeper_tb_2 = exception_tb;
    exception_keeper_lineno_2 = exception_lineno;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    exception_lineno = 0;

    // Preserve existing published exception id 1.
    GET_CURRENT_EXCEPTION(&exception_preserved_type_1, &exception_preserved_value_1, &exception_preserved_tb_1);

    if (exception_keeper_tb_2 == NULL) {
        exception_keeper_tb_2 = MAKE_TRACEBACK(frame_7fdb2c580e4cfeb4386220f78b54ec88, exception_keeper_lineno_2);
    } else if (exception_keeper_lineno_2 != 0) {
        exception_keeper_tb_2 = ADD_TRACEBACK(exception_keeper_tb_2, frame_7fdb2c580e4cfeb4386220f78b54ec88, exception_keeper_lineno_2);
    }

    NORMALIZE_EXCEPTION(&exception_keeper_type_2, &exception_keeper_value_2, &exception_keeper_tb_2);
    ATTACH_TRACEBACK_TO_EXCEPTION_VALUE(exception_keeper_value_2, exception_keeper_tb_2);
    PUBLISH_EXCEPTION(&exception_keeper_type_2, &exception_keeper_value_2, &exception_keeper_tb_2);
    // Tried code:
    {
        bool tmp_condition_result_2;
        PyObject *tmp_cmp_expr_left_1;
        PyObject *tmp_cmp_expr_right_1;
        tmp_cmp_expr_left_1 = EXC_TYPE(PyThreadState_GET());
        tmp_cmp_expr_right_1 = PyExc_BaseException;
        tmp_res = EXCEPTION_MATCH_BOOL(tmp_cmp_expr_left_1, tmp_cmp_expr_right_1);
        assert(!(tmp_res == -1));
        tmp_condition_result_2 = (tmp_res != 0) ? true : false;
        if (tmp_condition_result_2 != false) {
            goto branch_yes_2;
        } else {
            goto branch_no_2;
        }
    }
    branch_yes_2:;
    {
        nuitka_bool tmp_assign_source_31;
        tmp_assign_source_31 = NUITKA_BOOL_FALSE;
        tmp_with_1__indicator = tmp_assign_source_31;
    }
    {
        bool tmp_condition_result_3;
        PyObject *tmp_operand_value_1;
        PyObject *tmp_called_value_20;
        PyObject *tmp_args_element_value_11;
        PyObject *tmp_args_element_value_12;
        PyObject *tmp_args_element_value_13;
        CHECK_OBJECT(tmp_with_1__exit);
        tmp_called_value_20 = tmp_with_1__exit;
        tmp_args_element_value_11 = EXC_TYPE(PyThreadState_GET());
        tmp_args_element_value_12 = EXC_VALUE(PyThreadState_GET());
        tmp_args_element_value_13 = EXC_TRACEBACK(PyThreadState_GET());
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 66;
        {
            PyObject *call_args[] = {tmp_args_element_value_11, tmp_args_element_value_12, tmp_args_element_value_13};
            tmp_operand_value_1 = CALL_FUNCTION_WITH_ARGS3(tmp_called_value_20, call_args);
        }

        if (tmp_operand_value_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 66;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_7;
        }
        tmp_res = CHECK_IF_TRUE(tmp_operand_value_1);
        Py_DECREF(tmp_operand_value_1);
        if (tmp_res == -1) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 66;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_7;
        }
        tmp_condition_result_3 = (tmp_res == 0) ? true : false;
        if (tmp_condition_result_3 != false) {
            goto branch_yes_3;
        } else {
            goto branch_no_3;
        }
    }
    branch_yes_3:;
    tmp_result = RERAISE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
    if (unlikely(tmp_result == false)) {
        exception_lineno = 66;
    }

    if (exception_tb && exception_tb->tb_frame == &frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame) frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = exception_tb->tb_lineno;
    type_description_1 = "oooooooooooooooooooooo";
    goto try_except_handler_7;
    branch_no_3:;
    goto branch_end_2;
    branch_no_2:;
    tmp_result = RERAISE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
    if (unlikely(tmp_result == false)) {
        exception_lineno = 65;
    }

    if (exception_tb && exception_tb->tb_frame == &frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame) frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = exception_tb->tb_lineno;
    type_description_1 = "oooooooooooooooooooooo";
    goto try_except_handler_7;
    branch_end_2:;
    goto try_end_3;
    // Exception handler code:
    try_except_handler_7:;
    exception_keeper_type_3 = exception_type;
    exception_keeper_value_3 = exception_value;
    exception_keeper_tb_3 = exception_tb;
    exception_keeper_lineno_3 = exception_lineno;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    exception_lineno = 0;

    // Restore previous exception id 1.
    SET_CURRENT_EXCEPTION(exception_preserved_type_1, exception_preserved_value_1, exception_preserved_tb_1);

    // Re-raise.
    exception_type = exception_keeper_type_3;
    exception_value = exception_keeper_value_3;
    exception_tb = exception_keeper_tb_3;
    exception_lineno = exception_keeper_lineno_3;

    goto try_except_handler_4;
    // End of try:
    try_end_3:;
    // Restore previous exception id 1.
    SET_CURRENT_EXCEPTION(exception_preserved_type_1, exception_preserved_value_1, exception_preserved_tb_1);

    goto try_end_2;
    NUITKA_CANNOT_GET_HERE("exception handler codes exits in all cases");
    return NULL;
    // End of try:
    try_end_2:;
    goto try_end_4;
    // Exception handler code:
    try_except_handler_4:;
    exception_keeper_type_4 = exception_type;
    exception_keeper_value_4 = exception_value;
    exception_keeper_tb_4 = exception_tb;
    exception_keeper_lineno_4 = exception_lineno;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    exception_lineno = 0;

    {
        bool tmp_condition_result_4;
        nuitka_bool tmp_cmp_expr_left_2;
        nuitka_bool tmp_cmp_expr_right_2;
        assert(tmp_with_1__indicator != NUITKA_BOOL_UNASSIGNED);
        tmp_cmp_expr_left_2 = tmp_with_1__indicator;
        tmp_cmp_expr_right_2 = NUITKA_BOOL_TRUE;
        tmp_condition_result_4 = (tmp_cmp_expr_left_2 == tmp_cmp_expr_right_2) ? true : false;
        if (tmp_condition_result_4 != false) {
            goto branch_yes_4;
        } else {
            goto branch_no_4;
        }
    }
    branch_yes_4:;
    {
        PyObject *tmp_called_value_21;
        PyObject *tmp_call_result_11;
        CHECK_OBJECT(tmp_with_1__exit);
        tmp_called_value_21 = tmp_with_1__exit;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 66;
        tmp_call_result_11 = CALL_FUNCTION_WITH_POSARGS3(tmp_called_value_21, mod_consts[55]);

        if (tmp_call_result_11 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);

            Py_DECREF(exception_keeper_type_4);
            Py_XDECREF(exception_keeper_value_4);
            Py_XDECREF(exception_keeper_tb_4);

            exception_lineno = 66;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_3;
        }
        Py_DECREF(tmp_call_result_11);
    }
    branch_no_4:;
    // Re-raise.
    exception_type = exception_keeper_type_4;
    exception_value = exception_keeper_value_4;
    exception_tb = exception_keeper_tb_4;
    exception_lineno = exception_keeper_lineno_4;

    goto try_except_handler_3;
    // End of try:
    try_end_4:;
    {
        bool tmp_condition_result_5;
        nuitka_bool tmp_cmp_expr_left_3;
        nuitka_bool tmp_cmp_expr_right_3;
        assert(tmp_with_1__indicator != NUITKA_BOOL_UNASSIGNED);
        tmp_cmp_expr_left_3 = tmp_with_1__indicator;
        tmp_cmp_expr_right_3 = NUITKA_BOOL_TRUE;
        tmp_condition_result_5 = (tmp_cmp_expr_left_3 == tmp_cmp_expr_right_3) ? true : false;
        if (tmp_condition_result_5 != false) {
            goto branch_yes_5;
        } else {
            goto branch_no_5;
        }
    }
    branch_yes_5:;
    {
        PyObject *tmp_called_value_22;
        PyObject *tmp_call_result_12;
        CHECK_OBJECT(tmp_with_1__exit);
        tmp_called_value_22 = tmp_with_1__exit;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 66;
        tmp_call_result_12 = CALL_FUNCTION_WITH_POSARGS3(tmp_called_value_22, mod_consts[55]);

        if (tmp_call_result_12 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 66;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_3;
        }
        Py_DECREF(tmp_call_result_12);
    }
    branch_no_5:;
    goto try_end_5;
    // Exception handler code:
    try_except_handler_3:;
    exception_keeper_type_5 = exception_type;
    exception_keeper_value_5 = exception_value;
    exception_keeper_tb_5 = exception_tb;
    exception_keeper_lineno_5 = exception_lineno;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    exception_lineno = 0;

    Py_XDECREF(tmp_with_1__source);
    tmp_with_1__source = NULL;
    Py_XDECREF(tmp_with_1__enter);
    tmp_with_1__enter = NULL;
    Py_XDECREF(tmp_with_1__exit);
    tmp_with_1__exit = NULL;
    // Re-raise.
    exception_type = exception_keeper_type_5;
    exception_value = exception_keeper_value_5;
    exception_tb = exception_keeper_tb_5;
    exception_lineno = exception_keeper_lineno_5;

    goto try_except_handler_2;
    // End of try:
    try_end_5:;
    CHECK_OBJECT(tmp_with_1__source);
    Py_DECREF(tmp_with_1__source);
    tmp_with_1__source = NULL;
    CHECK_OBJECT(tmp_with_1__enter);
    Py_DECREF(tmp_with_1__enter);
    tmp_with_1__enter = NULL;
    Py_XDECREF(tmp_with_1__exit);
    tmp_with_1__exit = NULL;
    {
        PyObject *tmp_assign_source_32;
        tmp_assign_source_32 = PyList_New(0);
        {
            PyObject *old = var_unique_kmer_frequency;
            var_unique_kmer_frequency = tmp_assign_source_32;
            Py_XDECREF(old);
        }

    }
    {
        PyObject *tmp_called_value_23;
        PyObject *tmp_expression_value_25;
        PyObject *tmp_call_result_13;
        PyObject *tmp_args_element_value_14;
        CHECK_OBJECT(var_unique_kmer_frequency);
        tmp_expression_value_25 = var_unique_kmer_frequency;
        tmp_called_value_23 = LOOKUP_ATTRIBUTE(tmp_expression_value_25, mod_consts[56]);
        assert(!(tmp_called_value_23 == NULL));
        // Tried code:
        {
            PyObject *tmp_assign_source_33;
            PyObject *tmp_iter_arg_3;
            if (var_unique_sequences == NULL) {

                FORMAT_UNBOUND_CLOSURE_ERROR(&exception_type, &exception_value, mod_consts[57]);
                exception_tb = NULL;
                NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
                CHAIN_EXCEPTION(exception_value);

                exception_lineno = 71;
                type_description_1 = "oooooooooooooooooooooo";
                goto try_except_handler_8;
            }

            tmp_iter_arg_3 = var_unique_sequences;
            tmp_assign_source_33 = MAKE_ITERATOR(tmp_iter_arg_3);
            if (tmp_assign_source_33 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 71;
                type_description_1 = "oooooooooooooooooooooo";
                goto try_except_handler_8;
            }
            {
                PyObject *old = tmp_listcomp_1__$0;
                tmp_listcomp_1__$0 = tmp_assign_source_33;
                Py_XDECREF(old);
            }

        }
        {
            PyObject *tmp_assign_source_34;
            tmp_assign_source_34 = PyList_New(0);
            {
                PyObject *old = tmp_listcomp_1__contraction;
                tmp_listcomp_1__contraction = tmp_assign_source_34;
                Py_XDECREF(old);
            }

        }
        if (isFrameUnusable(cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2)) {
            Py_XDECREF(cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2);

#if _DEBUG_REFCOUNTS
            if (cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2 == NULL) {
                count_active_frame_cache_instances += 1;
            } else {
                count_released_frame_cache_instances += 1;
            }
            count_allocated_frame_cache_instances += 1;
#endif
            cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2 = MAKE_FUNCTION_FRAME(codeobj_ccd1fb807920dfdbbfd5c88e344983c6, module___main__, sizeof(void *)+sizeof(void *)+sizeof(void *));
#if _DEBUG_REFCOUNTS
        } else {
            count_hit_frame_cache_instances += 1;
#endif
        }
        assert(cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2->m_type_description == NULL);
        frame_ccd1fb807920dfdbbfd5c88e344983c6_2 = cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2;

        // Push the new frame as the currently active one.
        pushFrameStack(frame_ccd1fb807920dfdbbfd5c88e344983c6_2);

        // Mark the frame object as in use, ref count 1 will be up for reuse.
        assert(Py_REFCNT(frame_ccd1fb807920dfdbbfd5c88e344983c6_2) == 2); // Frame stack

        // Framed code:
        // Tried code:
        loop_start_3:;
        {
            PyObject *tmp_next_source_3;
            PyObject *tmp_assign_source_35;
            CHECK_OBJECT(tmp_listcomp_1__$0);
            tmp_next_source_3 = tmp_listcomp_1__$0;
            tmp_assign_source_35 = ITERATOR_NEXT(tmp_next_source_3);
            if (tmp_assign_source_35 == NULL) {
                if (CHECK_AND_CLEAR_STOP_ITERATION_OCCURRED()) {

                    goto loop_end_3;
                } else {

                    FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
                    type_description_2 = "ooo";
                    exception_lineno = 71;
                    goto try_except_handler_9;
                }
            }

            {
                PyObject *old = tmp_listcomp_1__iter_value_0;
                tmp_listcomp_1__iter_value_0 = tmp_assign_source_35;
                Py_XDECREF(old);
            }

        }
        {
            PyObject *tmp_assign_source_36;
            CHECK_OBJECT(tmp_listcomp_1__iter_value_0);
            tmp_assign_source_36 = tmp_listcomp_1__iter_value_0;
            {
                PyObject *old = outline_0_var_i;
                outline_0_var_i = tmp_assign_source_36;
                Py_INCREF(outline_0_var_i);
                Py_XDECREF(old);
            }

        }
        {
            PyObject *tmp_append_list_1;
            PyObject *tmp_append_value_1;
            PyObject *tmp_called_value_24;
            PyObject *tmp_expression_value_26;
            PyObject *tmp_expression_value_27;
            PyObject *tmp_subscript_value_7;
            PyObject *tmp_args_element_value_15;
            CHECK_OBJECT(tmp_listcomp_1__contraction);
            tmp_append_list_1 = tmp_listcomp_1__contraction;
            CHECK_OBJECT(var_input_sequence_list);
            tmp_expression_value_27 = var_input_sequence_list;
            CHECK_OBJECT(var_j);
            tmp_subscript_value_7 = var_j;
            tmp_expression_value_26 = LOOKUP_SUBSCRIPT(tmp_expression_value_27, tmp_subscript_value_7);
            if (tmp_expression_value_26 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 71;
                type_description_2 = "ooo";
                goto try_except_handler_9;
            }
            tmp_called_value_24 = LOOKUP_ATTRIBUTE(tmp_expression_value_26, mod_consts[58]);
            Py_DECREF(tmp_expression_value_26);
            if (tmp_called_value_24 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 71;
                type_description_2 = "ooo";
                goto try_except_handler_9;
            }
            CHECK_OBJECT(outline_0_var_i);
            tmp_args_element_value_15 = outline_0_var_i;
            frame_ccd1fb807920dfdbbfd5c88e344983c6_2->m_frame.f_lineno = 71;
            tmp_append_value_1 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_24, tmp_args_element_value_15);
            Py_DECREF(tmp_called_value_24);
            if (tmp_append_value_1 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 71;
                type_description_2 = "ooo";
                goto try_except_handler_9;
            }
            assert(PyList_Check(tmp_append_list_1));
            tmp_result = LIST_APPEND1(tmp_append_list_1, tmp_append_value_1);
            if (tmp_result == false) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 71;
                type_description_2 = "ooo";
                goto try_except_handler_9;
            }
        }
        if (CONSIDER_THREADING() == false) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 71;
            type_description_2 = "ooo";
            goto try_except_handler_9;
        }
        goto loop_start_3;
        loop_end_3:;
        CHECK_OBJECT(tmp_listcomp_1__contraction);
        tmp_args_element_value_14 = tmp_listcomp_1__contraction;
        Py_INCREF(tmp_args_element_value_14);
        goto try_return_handler_9;
        NUITKA_CANNOT_GET_HERE("tried codes exits in all cases");
        return NULL;
        // Return handler code:
        try_return_handler_9:;
        CHECK_OBJECT(tmp_listcomp_1__$0);
        Py_DECREF(tmp_listcomp_1__$0);
        tmp_listcomp_1__$0 = NULL;
        CHECK_OBJECT(tmp_listcomp_1__contraction);
        Py_DECREF(tmp_listcomp_1__contraction);
        tmp_listcomp_1__contraction = NULL;
        Py_XDECREF(tmp_listcomp_1__iter_value_0);
        tmp_listcomp_1__iter_value_0 = NULL;
        goto frame_return_exit_1;
        // Exception handler code:
        try_except_handler_9:;
        exception_keeper_type_6 = exception_type;
        exception_keeper_value_6 = exception_value;
        exception_keeper_tb_6 = exception_tb;
        exception_keeper_lineno_6 = exception_lineno;
        exception_type = NULL;
        exception_value = NULL;
        exception_tb = NULL;
        exception_lineno = 0;

        CHECK_OBJECT(tmp_listcomp_1__$0);
        Py_DECREF(tmp_listcomp_1__$0);
        tmp_listcomp_1__$0 = NULL;
        CHECK_OBJECT(tmp_listcomp_1__contraction);
        Py_DECREF(tmp_listcomp_1__contraction);
        tmp_listcomp_1__contraction = NULL;
        Py_XDECREF(tmp_listcomp_1__iter_value_0);
        tmp_listcomp_1__iter_value_0 = NULL;
        // Re-raise.
        exception_type = exception_keeper_type_6;
        exception_value = exception_keeper_value_6;
        exception_tb = exception_keeper_tb_6;
        exception_lineno = exception_keeper_lineno_6;

        goto frame_exception_exit_2;
        // End of try:

#if 0
        RESTORE_FRAME_EXCEPTION(frame_ccd1fb807920dfdbbfd5c88e344983c6_2);
#endif

        // Put the previous frame back on top.
        popFrameStack();

        goto frame_no_exception_1;

        frame_return_exit_1:;
#if 0
        RESTORE_FRAME_EXCEPTION(frame_ccd1fb807920dfdbbfd5c88e344983c6_2);
#endif

        // Put the previous frame back on top.
        popFrameStack();

        goto try_return_handler_8;

        frame_exception_exit_2:;

#if 0
        RESTORE_FRAME_EXCEPTION(frame_ccd1fb807920dfdbbfd5c88e344983c6_2);
#endif

        if (exception_tb == NULL) {
            exception_tb = MAKE_TRACEBACK(frame_ccd1fb807920dfdbbfd5c88e344983c6_2, exception_lineno);
        } else if (exception_tb->tb_frame != &frame_ccd1fb807920dfdbbfd5c88e344983c6_2->m_frame) {
            exception_tb = ADD_TRACEBACK(exception_tb, frame_ccd1fb807920dfdbbfd5c88e344983c6_2, exception_lineno);
        }

        // Attaches locals to frame if any.
        Nuitka_Frame_AttachLocals(
            frame_ccd1fb807920dfdbbfd5c88e344983c6_2,
            type_description_2,
            outline_0_var_i,
            var_input_sequence_list,
            var_j
        );


        // Release cached frame if used for exception.
        if (frame_ccd1fb807920dfdbbfd5c88e344983c6_2 == cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2) {
#if _DEBUG_REFCOUNTS
            count_active_frame_cache_instances -= 1;
            count_released_frame_cache_instances += 1;
#endif

            Py_DECREF(cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2);
            cache_frame_ccd1fb807920dfdbbfd5c88e344983c6_2 = NULL;
        }

        assertFrameObject(frame_ccd1fb807920dfdbbfd5c88e344983c6_2);

        // Put the previous frame back on top.
        popFrameStack();

        // Return the error.
        goto nested_frame_exit_1;

        frame_no_exception_1:;
        goto skip_nested_handling_1;
        nested_frame_exit_1:;
        type_description_1 = "oooooooooooooooooooooo";
        goto try_except_handler_8;
        skip_nested_handling_1:;
        NUITKA_CANNOT_GET_HERE("tried codes exits in all cases");
        return NULL;
        // Return handler code:
        try_return_handler_8:;
        Py_XDECREF(outline_0_var_i);
        outline_0_var_i = NULL;
        goto outline_result_1;
        // Exception handler code:
        try_except_handler_8:;
        exception_keeper_type_7 = exception_type;
        exception_keeper_value_7 = exception_value;
        exception_keeper_tb_7 = exception_tb;
        exception_keeper_lineno_7 = exception_lineno;
        exception_type = NULL;
        exception_value = NULL;
        exception_tb = NULL;
        exception_lineno = 0;

        Py_XDECREF(outline_0_var_i);
        outline_0_var_i = NULL;
        // Re-raise.
        exception_type = exception_keeper_type_7;
        exception_value = exception_keeper_value_7;
        exception_tb = exception_keeper_tb_7;
        exception_lineno = exception_keeper_lineno_7;

        goto outline_exception_1;
        // End of try:
        NUITKA_CANNOT_GET_HERE("Return statement must have exited already.");
        return NULL;
        outline_exception_1:;
        exception_lineno = 71;
        goto try_except_handler_2;
        outline_result_1:;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 71;
        tmp_call_result_13 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_23, tmp_args_element_value_14);
        Py_DECREF(tmp_called_value_23);
        Py_DECREF(tmp_args_element_value_14);
        if (tmp_call_result_13 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 71;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        Py_DECREF(tmp_call_result_13);
    }
    {
        PyObject *tmp_called_value_25;
        PyObject *tmp_expression_value_28;
        PyObject *tmp_call_result_14;
        PyObject *tmp_args_element_value_16;
        PyObject *tmp_called_value_26;
        PyObject *tmp_expression_value_29;
        PyObject *tmp_args_element_value_17;
        PyObject *tmp_args_element_value_18;
        if (var_prediction_output == NULL) {

            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[59]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 73;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_expression_value_28 = var_prediction_output;
        tmp_called_value_25 = LOOKUP_ATTRIBUTE(tmp_expression_value_28, mod_consts[56]);
        if (tmp_called_value_25 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 73;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        if (var_Prediction_object == NULL) {
            Py_DECREF(tmp_called_value_25);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[60]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 73;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_expression_value_29 = var_Prediction_object;
        tmp_called_value_26 = LOOKUP_ATTRIBUTE(tmp_expression_value_29, mod_consts[61]);
        if (tmp_called_value_26 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_25);

            exception_lineno = 73;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        CHECK_OBJECT(var_unique_kmer_frequency);
        tmp_args_element_value_17 = var_unique_kmer_frequency;
        if (var_plasmid == NULL) {
            Py_DECREF(tmp_called_value_25);
            Py_DECREF(tmp_called_value_26);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[19]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 73;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_args_element_value_18 = var_plasmid;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 73;
        {
            PyObject *call_args[] = {tmp_args_element_value_17, tmp_args_element_value_18};
            tmp_args_element_value_16 = CALL_FUNCTION_WITH_ARGS2(tmp_called_value_26, call_args);
        }

        Py_DECREF(tmp_called_value_26);
        if (tmp_args_element_value_16 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_25);

            exception_lineno = 73;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 73;
        tmp_call_result_14 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_25, tmp_args_element_value_16);
        Py_DECREF(tmp_called_value_25);
        Py_DECREF(tmp_args_element_value_16);
        if (tmp_call_result_14 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 73;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        Py_DECREF(tmp_call_result_14);
    }
    {
        PyObject *tmp_called_value_27;
        PyObject *tmp_expression_value_30;
        PyObject *tmp_call_result_15;
        PyObject *tmp_args_element_value_19;
        PyObject *tmp_called_value_28;
        PyObject *tmp_expression_value_31;
        PyObject *tmp_args_element_value_20;
        PyObject *tmp_called_value_29;
        PyObject *tmp_expression_value_32;
        PyObject *tmp_args_element_value_21;
        PyObject *tmp_args_element_value_22;
        PyObject *tmp_args_element_value_23;
        if (var_significance_output == NULL) {

            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[62]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_expression_value_30 = var_significance_output;
        tmp_called_value_27 = LOOKUP_ATTRIBUTE(tmp_expression_value_30, mod_consts[56]);
        if (tmp_called_value_27 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        if (var_Prediction_object == NULL) {
            Py_DECREF(tmp_called_value_27);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[60]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_expression_value_31 = var_Prediction_object;
        tmp_called_value_28 = LOOKUP_ATTRIBUTE(tmp_expression_value_31, mod_consts[63]);
        if (tmp_called_value_28 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_27);

            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        if (var_Prediction_object == NULL) {
            Py_DECREF(tmp_called_value_27);
            Py_DECREF(tmp_called_value_28);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[60]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_expression_value_32 = var_Prediction_object;
        tmp_called_value_29 = LOOKUP_ATTRIBUTE(tmp_expression_value_32, mod_consts[61]);
        if (tmp_called_value_29 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_27);
            Py_DECREF(tmp_called_value_28);

            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        CHECK_OBJECT(var_unique_kmer_frequency);
        tmp_args_element_value_21 = var_unique_kmer_frequency;
        if (var_plasmid == NULL) {
            Py_DECREF(tmp_called_value_27);
            Py_DECREF(tmp_called_value_28);
            Py_DECREF(tmp_called_value_29);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[19]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_args_element_value_22 = var_plasmid;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 74;
        {
            PyObject *call_args[] = {tmp_args_element_value_21, tmp_args_element_value_22};
            tmp_args_element_value_20 = CALL_FUNCTION_WITH_ARGS2(tmp_called_value_29, call_args);
        }

        Py_DECREF(tmp_called_value_29);
        if (tmp_args_element_value_20 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_27);
            Py_DECREF(tmp_called_value_28);

            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        if (var_plasmid == NULL) {
            Py_DECREF(tmp_called_value_27);
            Py_DECREF(tmp_called_value_28);
            Py_DECREF(tmp_args_element_value_20);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[19]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_args_element_value_23 = var_plasmid;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 74;
        {
            PyObject *call_args[] = {tmp_args_element_value_20, tmp_args_element_value_23};
            tmp_args_element_value_19 = CALL_FUNCTION_WITH_ARGS2(tmp_called_value_28, call_args);
        }

        Py_DECREF(tmp_called_value_28);
        Py_DECREF(tmp_args_element_value_20);
        if (tmp_args_element_value_19 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_27);

            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 74;
        tmp_call_result_15 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_27, tmp_args_element_value_19);
        Py_DECREF(tmp_called_value_27);
        Py_DECREF(tmp_args_element_value_19);
        if (tmp_call_result_15 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 74;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        Py_DECREF(tmp_call_result_15);
    }
    {
        nuitka_bool tmp_condition_result_6;
        int tmp_truth_name_2;
        CHECK_OBJECT(var_tree);
        tmp_truth_name_2 = CHECK_IF_TRUE(var_tree);
        if (tmp_truth_name_2 == -1) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 76;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        tmp_condition_result_6 = tmp_truth_name_2 == 0 ? NUITKA_BOOL_FALSE : NUITKA_BOOL_TRUE;
        if (tmp_condition_result_6 == NUITKA_BOOL_TRUE) {
            goto branch_yes_6;
        } else {
            goto branch_no_6;
        }
    }
    branch_yes_6:;
    {
        PyObject *tmp_assign_source_37;
        PyObject *tmp_left_value_5;
        PyObject *tmp_right_value_5;
        PyObject *tmp_called_value_30;
        PyObject *tmp_expression_value_33;
        PyObject *tmp_args_element_value_24;
        PyObject *tmp_expression_value_34;
        PyObject *tmp_subscript_value_8;
        PyObject *tmp_args_element_value_25;
        PyObject *tmp_expression_value_35;
        PyObject *tmp_subscript_value_9;
        if (var_converted_sequences_phyolgenetic_tree == NULL) {

            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[64]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 77;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_left_value_5 = var_converted_sequences_phyolgenetic_tree;
        if (var_Phylogenetic_object == NULL) {

            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[65]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 77;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_expression_value_33 = var_Phylogenetic_object;
        tmp_called_value_30 = LOOKUP_ATTRIBUTE(tmp_expression_value_33, mod_consts[66]);
        if (tmp_called_value_30 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 77;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        CHECK_OBJECT(var_unique_kmer_frequency);
        tmp_expression_value_34 = var_unique_kmer_frequency;
        tmp_subscript_value_8 = mod_consts[67];
        tmp_args_element_value_24 = LOOKUP_SUBSCRIPT_CONST(tmp_expression_value_34, tmp_subscript_value_8, 0);
        if (tmp_args_element_value_24 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_30);

            exception_lineno = 77;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        CHECK_OBJECT(var_label_list);
        tmp_expression_value_35 = var_label_list;
        CHECK_OBJECT(var_j);
        tmp_subscript_value_9 = var_j;
        tmp_args_element_value_25 = LOOKUP_SUBSCRIPT(tmp_expression_value_35, tmp_subscript_value_9);
        if (tmp_args_element_value_25 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_30);
            Py_DECREF(tmp_args_element_value_24);

            exception_lineno = 77;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 77;
        {
            PyObject *call_args[] = {tmp_args_element_value_24, tmp_args_element_value_25};
            tmp_right_value_5 = CALL_FUNCTION_WITH_ARGS2(tmp_called_value_30, call_args);
        }

        Py_DECREF(tmp_called_value_30);
        Py_DECREF(tmp_args_element_value_24);
        Py_DECREF(tmp_args_element_value_25);
        if (tmp_right_value_5 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 77;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        tmp_assign_source_37 = BINARY_OPERATION_ADD_OBJECT_OBJECT_OBJECT(tmp_left_value_5, tmp_right_value_5);
        Py_DECREF(tmp_right_value_5);
        if (tmp_assign_source_37 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 77;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        {
            PyObject *old = var_converted_sequences_phyolgenetic_tree;
            var_converted_sequences_phyolgenetic_tree = tmp_assign_source_37;
            Py_XDECREF(old);
        }

    }
    branch_no_6:;
    {
        PyObject *tmp_called_value_31;
        PyObject *tmp_expression_value_36;
        PyObject *tmp_call_result_16;
        PyObject *tmp_args_element_value_26;
        PyObject *tmp_called_value_32;
        PyObject *tmp_expression_value_37;
        PyObject *tmp_args_element_value_27;
        PyObject *tmp_args_element_value_28;
        PyObject *tmp_args_element_value_29;
        if (var_taxa_output == NULL) {

            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[68]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 78;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_expression_value_36 = var_taxa_output;
        tmp_called_value_31 = LOOKUP_ATTRIBUTE(tmp_expression_value_36, mod_consts[56]);
        if (tmp_called_value_31 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 78;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        if (var_Prediction_object == NULL) {
            Py_DECREF(tmp_called_value_31);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[60]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 78;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_expression_value_37 = var_Prediction_object;
        tmp_called_value_32 = LOOKUP_ATTRIBUTE(tmp_expression_value_37, mod_consts[69]);
        if (tmp_called_value_32 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_31);

            exception_lineno = 78;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        CHECK_OBJECT(var_unique_kmer_frequency);
        tmp_args_element_value_27 = var_unique_kmer_frequency;
        if (var_plasmid == NULL) {
            Py_DECREF(tmp_called_value_31);
            Py_DECREF(tmp_called_value_32);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[19]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 78;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_args_element_value_28 = var_plasmid;
        if (var_rank == NULL) {
            Py_DECREF(tmp_called_value_31);
            Py_DECREF(tmp_called_value_32);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[20]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 78;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }

        tmp_args_element_value_29 = var_rank;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 78;
        {
            PyObject *call_args[] = {tmp_args_element_value_27, tmp_args_element_value_28, tmp_args_element_value_29};
            tmp_args_element_value_26 = CALL_FUNCTION_WITH_ARGS3(tmp_called_value_32, call_args);
        }

        Py_DECREF(tmp_called_value_32);
        if (tmp_args_element_value_26 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_31);

            exception_lineno = 78;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 78;
        tmp_call_result_16 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_31, tmp_args_element_value_26);
        Py_DECREF(tmp_called_value_31);
        Py_DECREF(tmp_args_element_value_26);
        if (tmp_call_result_16 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 78;
            type_description_1 = "oooooooooooooooooooooo";
            goto try_except_handler_2;
        }
        Py_DECREF(tmp_call_result_16);
    }
    if (CONSIDER_THREADING() == false) {
        assert(ERROR_OCCURRED());

        FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


        exception_lineno = 63;
        type_description_1 = "oooooooooooooooooooooo";
        goto try_except_handler_2;
    }
    goto loop_start_1;
    loop_end_1:;
    goto try_end_6;
    // Exception handler code:
    try_except_handler_2:;
    exception_keeper_type_8 = exception_type;
    exception_keeper_value_8 = exception_value;
    exception_keeper_tb_8 = exception_tb;
    exception_keeper_lineno_8 = exception_lineno;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    exception_lineno = 0;

    Py_XDECREF(tmp_for_loop_1__iter_value);
    tmp_for_loop_1__iter_value = NULL;
    CHECK_OBJECT(tmp_for_loop_1__for_iterator);
    Py_DECREF(tmp_for_loop_1__for_iterator);
    tmp_for_loop_1__for_iterator = NULL;
    // Re-raise.
    exception_type = exception_keeper_type_8;
    exception_value = exception_keeper_value_8;
    exception_tb = exception_keeper_tb_8;
    exception_lineno = exception_keeper_lineno_8;

    goto frame_exception_exit_1;
    // End of try:
    try_end_6:;
    Py_XDECREF(tmp_for_loop_1__iter_value);
    tmp_for_loop_1__iter_value = NULL;
    CHECK_OBJECT(tmp_for_loop_1__for_iterator);
    Py_DECREF(tmp_for_loop_1__for_iterator);
    tmp_for_loop_1__for_iterator = NULL;
    {
        PyObject *tmp_assign_source_38;
        PyObject *tmp_called_value_33;
        PyObject *tmp_expression_value_38;
        PyObject *tmp_args_value_3;
        PyObject *tmp_tuple_element_3;
        PyObject *tmp_list_element_1;
        PyObject *tmp_called_instance_5;
        PyObject *tmp_args_element_value_30;
        PyObject *tmp_kwargs_value_3;
        tmp_expression_value_38 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[33]);

        if (unlikely(tmp_expression_value_38 == NULL)) {
            tmp_expression_value_38 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[33]);
        }

        if (tmp_expression_value_38 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 80;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_called_value_33 = LOOKUP_ATTRIBUTE(tmp_expression_value_38, mod_consts[70]);
        if (tmp_called_value_33 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 80;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_called_instance_5 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[33]);

        if (unlikely(tmp_called_instance_5 == NULL)) {
            tmp_called_instance_5 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[33]);
        }

        if (tmp_called_instance_5 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_33);

            exception_lineno = 80;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        CHECK_OBJECT(var_input_sequence_list);
        tmp_args_element_value_30 = var_input_sequence_list;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 80;
        tmp_list_element_1 = CALL_METHOD_WITH_SINGLE_ARG(tmp_called_instance_5, mod_consts[71], tmp_args_element_value_30);
        if (tmp_list_element_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_33);

            exception_lineno = 80;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_tuple_element_3 = PyList_New(4);
        {
            PyObject *tmp_called_value_34;
            PyObject *tmp_expression_value_39;
            PyObject *tmp_args_element_value_31;
            PyObject *tmp_called_value_35;
            PyObject *tmp_expression_value_40;
            PyObject *tmp_args_element_value_32;
            PyObject *tmp_called_value_36;
            PyObject *tmp_expression_value_41;
            PyObject *tmp_args_element_value_33;
            PyList_SET_ITEM(tmp_tuple_element_3, 0, tmp_list_element_1);
            tmp_expression_value_39 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[33]);

            if (unlikely(tmp_expression_value_39 == NULL)) {
                tmp_expression_value_39 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[33]);
            }

            if (tmp_expression_value_39 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }
            tmp_called_value_34 = LOOKUP_ATTRIBUTE(tmp_expression_value_39, mod_consts[71]);
            if (tmp_called_value_34 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }
            if (var_prediction_output == NULL) {
                Py_DECREF(tmp_called_value_34);
                FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[59]);
                exception_tb = NULL;
                NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
                CHAIN_EXCEPTION(exception_value);

                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }

            tmp_args_element_value_31 = var_prediction_output;
            frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 80;
            tmp_list_element_1 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_34, tmp_args_element_value_31);
            Py_DECREF(tmp_called_value_34);
            if (tmp_list_element_1 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }
            PyList_SET_ITEM(tmp_tuple_element_3, 1, tmp_list_element_1);
            tmp_expression_value_40 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[33]);

            if (unlikely(tmp_expression_value_40 == NULL)) {
                tmp_expression_value_40 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[33]);
            }

            if (tmp_expression_value_40 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }
            tmp_called_value_35 = LOOKUP_ATTRIBUTE(tmp_expression_value_40, mod_consts[71]);
            if (tmp_called_value_35 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }
            if (var_taxa_output == NULL) {
                Py_DECREF(tmp_called_value_35);
                FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[68]);
                exception_tb = NULL;
                NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
                CHAIN_EXCEPTION(exception_value);

                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }

            tmp_args_element_value_32 = var_taxa_output;
            frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 80;
            tmp_list_element_1 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_35, tmp_args_element_value_32);
            Py_DECREF(tmp_called_value_35);
            if (tmp_list_element_1 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }
            PyList_SET_ITEM(tmp_tuple_element_3, 2, tmp_list_element_1);
            tmp_expression_value_41 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[33]);

            if (unlikely(tmp_expression_value_41 == NULL)) {
                tmp_expression_value_41 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[33]);
            }

            if (tmp_expression_value_41 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }
            tmp_called_value_36 = LOOKUP_ATTRIBUTE(tmp_expression_value_41, mod_consts[71]);
            if (tmp_called_value_36 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }
            if (var_significance_output == NULL) {
                Py_DECREF(tmp_called_value_36);
                FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[62]);
                exception_tb = NULL;
                NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
                CHAIN_EXCEPTION(exception_value);

                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }

            tmp_args_element_value_33 = var_significance_output;
            frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 80;
            tmp_list_element_1 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_36, tmp_args_element_value_33);
            Py_DECREF(tmp_called_value_36);
            if (tmp_list_element_1 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 80;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_1;
            }
            PyList_SET_ITEM(tmp_tuple_element_3, 3, tmp_list_element_1);
        }
        goto list_build_noexception_1;
        // Exception handling pass through code for list_build:
        list_build_exception_1:;
        Py_DECREF(tmp_called_value_33);
        Py_DECREF(tmp_tuple_element_3);
        goto frame_exception_exit_1;
        // Finished with no exception for list_build:
        list_build_noexception_1:;
        tmp_args_value_3 = PyTuple_New(1);
        PyTuple_SET_ITEM(tmp_args_value_3, 0, tmp_tuple_element_3);
        tmp_kwargs_value_3 = PyDict_Copy(mod_consts[72]);
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 80;
        tmp_assign_source_38 = CALL_FUNCTION(tmp_called_value_33, tmp_args_value_3, tmp_kwargs_value_3);
        Py_DECREF(tmp_called_value_33);
        Py_DECREF(tmp_args_value_3);
        Py_DECREF(tmp_kwargs_value_3);
        if (tmp_assign_source_38 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 80;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        {
            PyObject *old = var_prediction_output;
            var_prediction_output = tmp_assign_source_38;
            Py_XDECREF(old);
        }

    }
    {
        PyObject *tmp_assattr_value_1;
        PyObject *tmp_list_element_2;
        PyObject *tmp_assattr_target_1;
        tmp_list_element_2 = mod_consts[73];
        tmp_assattr_value_1 = PyList_New(4);
        {
            PyObject *tmp_left_value_6;
            PyObject *tmp_right_value_6;
            PyList_SET_ITEM0(tmp_assattr_value_1, 0, tmp_list_element_2);
            tmp_list_element_2 = mod_consts[19];
            PyList_SET_ITEM0(tmp_assattr_value_1, 1, tmp_list_element_2);
            tmp_left_value_6 = mod_consts[74];
            if (var_rank == NULL) {

                FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[20]);
                exception_tb = NULL;
                NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
                CHAIN_EXCEPTION(exception_value);

                exception_lineno = 81;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_2;
            }

            tmp_right_value_6 = var_rank;
            tmp_list_element_2 = BINARY_OPERATION_ADD_OBJECT_UNICODE_OBJECT(tmp_left_value_6, tmp_right_value_6);
            if (tmp_list_element_2 == NULL) {
                assert(ERROR_OCCURRED());

                FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


                exception_lineno = 81;
                type_description_1 = "oooooooooooooooooooooo";
                goto list_build_exception_2;
            }
            PyList_SET_ITEM(tmp_assattr_value_1, 2, tmp_list_element_2);
            tmp_list_element_2 = mod_consts[75];
            PyList_SET_ITEM0(tmp_assattr_value_1, 3, tmp_list_element_2);
        }
        goto list_build_noexception_2;
        // Exception handling pass through code for list_build:
        list_build_exception_2:;
        Py_DECREF(tmp_assattr_value_1);
        goto frame_exception_exit_1;
        // Finished with no exception for list_build:
        list_build_noexception_2:;
        CHECK_OBJECT(var_prediction_output);
        tmp_assattr_target_1 = var_prediction_output;
        tmp_result = SET_ATTRIBUTE(tmp_assattr_target_1, mod_consts[76], tmp_assattr_value_1);
        Py_DECREF(tmp_assattr_value_1);
        if (tmp_result == false) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 81;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
    }
    {
        PyObject *tmp_assattr_value_2;
        PyObject *tmp_assattr_target_2;
        if (var_plasmid == NULL) {

            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[19]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 82;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }

        tmp_assattr_value_2 = var_plasmid;
        CHECK_OBJECT(var_prediction_output);
        tmp_assattr_target_2 = var_prediction_output;
        tmp_result = SET_ATTRIBUTE(tmp_assattr_target_2, mod_consts[19], tmp_assattr_value_2);
        if (tmp_result == false) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 82;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
    }
    {
        PyObject *tmp_assattr_value_3;
        PyObject *tmp_assattr_target_3;
        CHECK_OBJECT(var_label_list);
        tmp_assattr_value_3 = var_label_list;
        CHECK_OBJECT(var_prediction_output);
        tmp_assattr_target_3 = var_prediction_output;
        tmp_result = SET_ATTRIBUTE(tmp_assattr_target_3, mod_consts[77], tmp_assattr_value_3);
        if (tmp_result == false) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 83;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
    }
    {
        PyObject *tmp_assign_source_39;
        PyObject *tmp_called_instance_6;
        CHECK_OBJECT(var_prediction_output);
        tmp_called_instance_6 = var_prediction_output;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 84;
        tmp_assign_source_39 = CALL_METHOD_WITH_SINGLE_ARG(
            tmp_called_instance_6,
            mod_consts[78],
            PyTuple_GET_ITEM(mod_consts[79], 0)
        );

        if (tmp_assign_source_39 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 84;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        {
            PyObject *old = var_prediction_output;
            assert(old != NULL);
            var_prediction_output = tmp_assign_source_39;
            Py_DECREF(old);
        }

    }
    {
        PyObject *tmp_called_value_37;
        PyObject *tmp_expression_value_42;
        PyObject *tmp_call_result_17;
        PyObject *tmp_args_element_value_34;
        PyObject *tmp_called_value_38;
        PyObject *tmp_expression_value_43;
        PyObject *tmp_expression_value_44;
        PyObject *tmp_args_element_value_35;
        PyObject *tmp_args_element_value_36;
        PyObject *tmp_left_value_7;
        PyObject *tmp_expression_value_45;
        PyObject *tmp_subscript_value_10;
        PyObject *tmp_right_value_7;
        CHECK_OBJECT(var_prediction_output);
        tmp_expression_value_42 = var_prediction_output;
        tmp_called_value_37 = LOOKUP_ATTRIBUTE(tmp_expression_value_42, mod_consts[80]);
        if (tmp_called_value_37 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 85;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_expression_value_44 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[27]);

        if (unlikely(tmp_expression_value_44 == NULL)) {
            tmp_expression_value_44 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[27]);
        }

        if (tmp_expression_value_44 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_37);

            exception_lineno = 85;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_expression_value_43 = LOOKUP_ATTRIBUTE(tmp_expression_value_44, mod_consts[28]);
        if (tmp_expression_value_43 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_37);

            exception_lineno = 85;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_called_value_38 = LOOKUP_ATTRIBUTE(tmp_expression_value_43, mod_consts[29]);
        Py_DECREF(tmp_expression_value_43);
        if (tmp_called_value_38 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_37);

            exception_lineno = 85;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_args_element_value_35 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[23]);

        if (unlikely(tmp_args_element_value_35 == NULL)) {
            tmp_args_element_value_35 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[23]);
        }

        if (tmp_args_element_value_35 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_37);
            Py_DECREF(tmp_called_value_38);

            exception_lineno = 85;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        CHECK_OBJECT(var_args);
        tmp_expression_value_45 = var_args;
        tmp_subscript_value_10 = mod_consts[81];
        tmp_left_value_7 = LOOKUP_SUBSCRIPT(tmp_expression_value_45, tmp_subscript_value_10);
        if (tmp_left_value_7 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_37);
            Py_DECREF(tmp_called_value_38);

            exception_lineno = 85;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_right_value_7 = mod_consts[37];
        tmp_args_element_value_36 = BINARY_OPERATION_ADD_OBJECT_OBJECT_UNICODE(tmp_left_value_7, tmp_right_value_7);
        Py_DECREF(tmp_left_value_7);
        if (tmp_args_element_value_36 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_37);
            Py_DECREF(tmp_called_value_38);

            exception_lineno = 85;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 85;
        {
            PyObject *call_args[] = {tmp_args_element_value_35, tmp_args_element_value_36};
            tmp_args_element_value_34 = CALL_FUNCTION_WITH_ARGS2(tmp_called_value_38, call_args);
        }

        Py_DECREF(tmp_called_value_38);
        Py_DECREF(tmp_args_element_value_36);
        if (tmp_args_element_value_34 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_37);

            exception_lineno = 85;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 85;
        tmp_call_result_17 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_37, tmp_args_element_value_34);
        Py_DECREF(tmp_called_value_37);
        Py_DECREF(tmp_args_element_value_34);
        if (tmp_call_result_17 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 85;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_17);
    }
    {
        nuitka_bool tmp_condition_result_7;
        int tmp_truth_name_3;
        CHECK_OBJECT(var_tree);
        tmp_truth_name_3 = CHECK_IF_TRUE(var_tree);
        if (tmp_truth_name_3 == -1) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 88;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        tmp_condition_result_7 = tmp_truth_name_3 == 0 ? NUITKA_BOOL_FALSE : NUITKA_BOOL_TRUE;
        if (tmp_condition_result_7 == NUITKA_BOOL_TRUE) {
            goto branch_yes_7;
        } else {
            goto branch_no_7;
        }
    }
    branch_yes_7:;
    {
        PyObject *tmp_called_value_39;
        PyObject *tmp_call_result_18;
        tmp_called_value_39 = LOOKUP_BUILTIN(mod_consts[25]);
        assert(tmp_called_value_39 != NULL);
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 89;
        tmp_call_result_18 = CALL_FUNCTION_WITH_POSARGS1(tmp_called_value_39, mod_consts[82]);

        if (tmp_call_result_18 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 89;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_18);
    }
    {
        PyObject *tmp_called_value_40;
        PyObject *tmp_expression_value_46;
        PyObject *tmp_call_result_19;
        PyObject *tmp_args_element_value_37;
        if (var_Phylogenetic_object == NULL) {

            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[65]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 90;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }

        tmp_expression_value_46 = var_Phylogenetic_object;
        tmp_called_value_40 = LOOKUP_ATTRIBUTE(tmp_expression_value_46, mod_consts[83]);
        if (tmp_called_value_40 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 90;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        if (var_converted_sequences_phyolgenetic_tree == NULL) {
            Py_DECREF(tmp_called_value_40);
            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[64]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 90;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }

        tmp_args_element_value_37 = var_converted_sequences_phyolgenetic_tree;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 90;
        tmp_call_result_19 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_40, tmp_args_element_value_37);
        Py_DECREF(tmp_called_value_40);
        if (tmp_call_result_19 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 90;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_19);
    }
    {
        PyObject *tmp_assign_source_40;
        PyObject *tmp_called_instance_7;
        if (var_Phylogenetic_object == NULL) {

            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[65]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 91;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }

        tmp_called_instance_7 = var_Phylogenetic_object;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 91;
        tmp_assign_source_40 = CALL_METHOD_NO_ARGS(tmp_called_instance_7, mod_consts[84]);
        if (tmp_assign_source_40 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 91;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        assert(var_distance_matrix == NULL);
        var_distance_matrix = tmp_assign_source_40;
    }
    {
        PyObject *tmp_called_instance_8;
        PyObject *tmp_call_result_20;
        PyObject *tmp_args_element_value_38;
        if (var_Phylogenetic_object == NULL) {

            FORMAT_UNBOUND_LOCAL_ERROR(&exception_type, &exception_value, mod_consts[65]);
            exception_tb = NULL;
            NORMALIZE_EXCEPTION(&exception_type, &exception_value, &exception_tb);
            CHAIN_EXCEPTION(exception_value);

            exception_lineno = 92;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }

        tmp_called_instance_8 = var_Phylogenetic_object;
        CHECK_OBJECT(var_distance_matrix);
        tmp_args_element_value_38 = var_distance_matrix;
        frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame.f_lineno = 92;
        tmp_call_result_20 = CALL_METHOD_WITH_SINGLE_ARG(tmp_called_instance_8, mod_consts[85], tmp_args_element_value_38);
        if (tmp_call_result_20 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 92;
            type_description_1 = "oooooooooooooooooooooo";
            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_20);
    }
    branch_no_7:;

#if 0
    RESTORE_FRAME_EXCEPTION(frame_7fdb2c580e4cfeb4386220f78b54ec88);
#endif

    // Put the previous frame back on top.
    popFrameStack();

    goto frame_no_exception_2;

    frame_exception_exit_1:;

#if 0
    RESTORE_FRAME_EXCEPTION(frame_7fdb2c580e4cfeb4386220f78b54ec88);
#endif

    if (exception_tb == NULL) {
        exception_tb = MAKE_TRACEBACK(frame_7fdb2c580e4cfeb4386220f78b54ec88, exception_lineno);
    } else if (exception_tb->tb_frame != &frame_7fdb2c580e4cfeb4386220f78b54ec88->m_frame) {
        exception_tb = ADD_TRACEBACK(exception_tb, frame_7fdb2c580e4cfeb4386220f78b54ec88, exception_lineno);
    }

    // Attaches locals to frame if any.
    Nuitka_Frame_AttachLocals(
        frame_7fdb2c580e4cfeb4386220f78b54ec88,
        type_description_1,
        var_parser,
        var_requiredNamed,
        var_args,
        var_plasmid,
        var_rank,
        var_tree,
        var_Prediction_object,
        var_Phylogenetic_object,
        var_input_file,
        var_files_lines,
        var_unique_sequences,
        var_input_sequence_list,
        var_label_list,
        var_prediction_output,
        var_significance_output,
        var_taxa_output,
        var_converted_sequences_phyolgenetic_tree,
        var_j,
        var_bar,
        var_user,
        var_unique_kmer_frequency,
        var_distance_matrix
    );


    // Release cached frame if used for exception.
    if (frame_7fdb2c580e4cfeb4386220f78b54ec88 == cache_frame_7fdb2c580e4cfeb4386220f78b54ec88) {
#if _DEBUG_REFCOUNTS
        count_active_frame_cache_instances -= 1;
        count_released_frame_cache_instances += 1;
#endif

        Py_DECREF(cache_frame_7fdb2c580e4cfeb4386220f78b54ec88);
        cache_frame_7fdb2c580e4cfeb4386220f78b54ec88 = NULL;
    }

    assertFrameObject(frame_7fdb2c580e4cfeb4386220f78b54ec88);

    // Put the previous frame back on top.
    popFrameStack();

    // Return the error.
    goto try_except_handler_1;

    frame_no_exception_2:;
    tmp_return_value = Py_None;
    Py_INCREF(tmp_return_value);
    goto try_return_handler_1;
    NUITKA_CANNOT_GET_HERE("tried codes exits in all cases");
    return NULL;
    // Return handler code:
    try_return_handler_1:;
    CHECK_OBJECT(var_parser);
    Py_DECREF(var_parser);
    var_parser = NULL;
    CHECK_OBJECT(var_requiredNamed);
    Py_DECREF(var_requiredNamed);
    var_requiredNamed = NULL;
    CHECK_OBJECT(var_args);
    Py_DECREF(var_args);
    var_args = NULL;
    Py_XDECREF(var_plasmid);
    var_plasmid = NULL;
    Py_XDECREF(var_rank);
    var_rank = NULL;
    CHECK_OBJECT(var_tree);
    Py_DECREF(var_tree);
    var_tree = NULL;
    Py_XDECREF(var_Prediction_object);
    var_Prediction_object = NULL;
    Py_XDECREF(var_Phylogenetic_object);
    var_Phylogenetic_object = NULL;
    CHECK_OBJECT(var_input_file);
    Py_DECREF(var_input_file);
    var_input_file = NULL;
    CHECK_OBJECT(var_files_lines);
    Py_DECREF(var_files_lines);
    var_files_lines = NULL;
    Py_XDECREF(var_unique_sequences);
    var_unique_sequences = NULL;
    CHECK_OBJECT(var_input_sequence_list);
    Py_DECREF(var_input_sequence_list);
    var_input_sequence_list = NULL;
    CHECK_OBJECT(var_label_list);
    Py_DECREF(var_label_list);
    var_label_list = NULL;
    CHECK_OBJECT(var_prediction_output);
    Py_DECREF(var_prediction_output);
    var_prediction_output = NULL;
    Py_XDECREF(var_significance_output);
    var_significance_output = NULL;
    Py_XDECREF(var_taxa_output);
    var_taxa_output = NULL;
    Py_XDECREF(var_converted_sequences_phyolgenetic_tree);
    var_converted_sequences_phyolgenetic_tree = NULL;
    Py_XDECREF(var_j);
    var_j = NULL;
    Py_XDECREF(var_bar);
    var_bar = NULL;
    Py_XDECREF(var_user);
    var_user = NULL;
    Py_XDECREF(var_unique_kmer_frequency);
    var_unique_kmer_frequency = NULL;
    Py_XDECREF(var_distance_matrix);
    var_distance_matrix = NULL;
    goto function_return_exit;
    // Exception handler code:
    try_except_handler_1:;
    exception_keeper_type_9 = exception_type;
    exception_keeper_value_9 = exception_value;
    exception_keeper_tb_9 = exception_tb;
    exception_keeper_lineno_9 = exception_lineno;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    exception_lineno = 0;

    Py_XDECREF(var_parser);
    var_parser = NULL;
    Py_XDECREF(var_requiredNamed);
    var_requiredNamed = NULL;
    Py_XDECREF(var_args);
    var_args = NULL;
    Py_XDECREF(var_plasmid);
    var_plasmid = NULL;
    Py_XDECREF(var_rank);
    var_rank = NULL;
    Py_XDECREF(var_tree);
    var_tree = NULL;
    Py_XDECREF(var_Prediction_object);
    var_Prediction_object = NULL;
    Py_XDECREF(var_Phylogenetic_object);
    var_Phylogenetic_object = NULL;
    Py_XDECREF(var_input_file);
    var_input_file = NULL;
    Py_XDECREF(var_files_lines);
    var_files_lines = NULL;
    Py_XDECREF(var_unique_sequences);
    var_unique_sequences = NULL;
    Py_XDECREF(var_input_sequence_list);
    var_input_sequence_list = NULL;
    Py_XDECREF(var_label_list);
    var_label_list = NULL;
    Py_XDECREF(var_prediction_output);
    var_prediction_output = NULL;
    Py_XDECREF(var_significance_output);
    var_significance_output = NULL;
    Py_XDECREF(var_taxa_output);
    var_taxa_output = NULL;
    Py_XDECREF(var_converted_sequences_phyolgenetic_tree);
    var_converted_sequences_phyolgenetic_tree = NULL;
    Py_XDECREF(var_j);
    var_j = NULL;
    Py_XDECREF(var_bar);
    var_bar = NULL;
    Py_XDECREF(var_user);
    var_user = NULL;
    Py_XDECREF(var_unique_kmer_frequency);
    var_unique_kmer_frequency = NULL;
    Py_XDECREF(var_distance_matrix);
    var_distance_matrix = NULL;
    // Re-raise.
    exception_type = exception_keeper_type_9;
    exception_value = exception_keeper_value_9;
    exception_tb = exception_keeper_tb_9;
    exception_lineno = exception_keeper_lineno_9;

    goto function_exception_exit;
    // End of try:

    NUITKA_CANNOT_GET_HERE("Return statement must have exited already.");
    return NULL;

function_exception_exit:
    assert(exception_type);
    RESTORE_ERROR_OCCURRED(exception_type, exception_value, exception_tb);

    return NULL;

function_return_exit:
   // Function cleanup code if any.


   // Actual function exit with return value, making sure we did not make
   // the error status worse despite non-NULL return.
   CHECK_OBJECT(tmp_return_value);
   assert(had_error || !ERROR_OCCURRED());
   return tmp_return_value;
}



static PyObject *MAKE_FUNCTION___main__$$$function__1_main() {
    struct Nuitka_FunctionObject *result = Nuitka_Function_New(
        impl___main__$$$function__1_main,
        mod_consts[104],
#if PYTHON_VERSION >= 0x300
        NULL,
#endif
        codeobj_7fdb2c580e4cfeb4386220f78b54ec88,
        NULL,
#if PYTHON_VERSION >= 0x300
        NULL,
        NULL,
#endif
        module___main__,
        NULL,
        NULL,
        0
    );


    return (PyObject *)result;
}


extern void _initCompiledCellType();
extern void _initCompiledGeneratorType();
extern void _initCompiledFunctionType();
extern void _initCompiledMethodType();
extern void _initCompiledFrameType();

extern PyTypeObject Nuitka_Loader_Type;

#ifdef _NUITKA_PLUGIN_DILL_ENABLED
// Provide a way to create find a function via its C code and create it back
// in another process, useful for multiprocessing extensions like dill
extern void registerDillPluginTables(char const *module_name, PyMethodDef *reduce_compiled_function, PyMethodDef *create_compiled_function);

function_impl_code functable___main__[] = {
    impl___main__$$$function__1_main,
    NULL
};

static char const *_reduce_compiled_function_argnames[] = {
    "func",
    NULL
};

static PyObject *_reduce_compiled_function(PyObject *self, PyObject *args, PyObject *kwds) {
    PyObject *func;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O:reduce_compiled_function", (char **)_reduce_compiled_function_argnames, &func, NULL)) {
        return NULL;
    }

    if (Nuitka_Function_Check(func) == false) {
        SET_CURRENT_EXCEPTION_TYPE0_STR(PyExc_TypeError, "not a compiled function");
        return NULL;
    }

    struct Nuitka_FunctionObject *function = (struct Nuitka_FunctionObject *)func;

    function_impl_code *current = functable___main__;
    int offset = 0;

    while (*current != NULL) {
        if (*current == function->m_c_code) {
            break;
        }

        current += 1;
        offset += 1;
    }

    if (*current == NULL) {
        SET_CURRENT_EXCEPTION_TYPE0_STR(PyExc_TypeError, "Cannot find compiled function in module.");
        return NULL;
    }

    PyObject *code_object_desc = PyTuple_New(6);
    PyTuple_SET_ITEM0(code_object_desc, 0, function->m_code_object->co_filename);
    PyTuple_SET_ITEM0(code_object_desc, 1, function->m_code_object->co_name);
    PyTuple_SET_ITEM(code_object_desc, 2, PyLong_FromLong(function->m_code_object->co_firstlineno));
    PyTuple_SET_ITEM0(code_object_desc, 3, function->m_code_object->co_varnames);
    PyTuple_SET_ITEM(code_object_desc, 4, PyLong_FromLong(function->m_code_object->co_argcount));
    PyTuple_SET_ITEM(code_object_desc, 5, PyLong_FromLong(function->m_code_object->co_flags));

    CHECK_OBJECT_DEEP(code_object_desc);

    PyObject *result = PyTuple_New(4);
    PyTuple_SET_ITEM(result, 0, PyLong_FromLong(offset));
    PyTuple_SET_ITEM(result, 1, code_object_desc);
    PyTuple_SET_ITEM0(result, 2, function->m_defaults);
    PyTuple_SET_ITEM0(result, 3, function->m_doc != NULL ? function->m_doc : Py_None);

    CHECK_OBJECT_DEEP(result);

    return result;
}

static PyMethodDef _method_def_reduce_compiled_function = {"reduce_compiled_function", (PyCFunction)_reduce_compiled_function,
                                                           METH_VARARGS | METH_KEYWORDS, NULL};

static char const *_create_compiled_function_argnames[] = {
    "func",
    "code_object_desc",
    "defaults",
    "doc",
    NULL
};


static PyObject *_create_compiled_function(PyObject *self, PyObject *args, PyObject *kwds) {
    CHECK_OBJECT_DEEP(args);

    PyObject *func;
    PyObject *code_object_desc;
    PyObject *defaults;
    PyObject *doc;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOO:create_compiled_function", (char **)_create_compiled_function_argnames, &func, &code_object_desc, &defaults, &doc, NULL)) {
        return NULL;
    }

    int offset = PyLong_AsLong(func);

    if (offset == -1 && ERROR_OCCURRED()) {
        return NULL;
    }

    if (offset > sizeof(functable___main__) || offset < 0) {
        SET_CURRENT_EXCEPTION_TYPE0_STR(PyExc_TypeError, "Wrong offset for compiled function.");
        return NULL;
    }

    PyObject *filename = PyTuple_GET_ITEM(code_object_desc, 0);
    PyObject *function_name = PyTuple_GET_ITEM(code_object_desc, 1);
    PyObject *line = PyTuple_GET_ITEM(code_object_desc, 2);
    int line_int = PyLong_AsLong(line);
    assert(!ERROR_OCCURRED());

    PyObject *argnames = PyTuple_GET_ITEM(code_object_desc, 3);
    PyObject *arg_count = PyTuple_GET_ITEM(code_object_desc, 4);
    int arg_count_int = PyLong_AsLong(arg_count);
    assert(!ERROR_OCCURRED());
    PyObject *flags = PyTuple_GET_ITEM(code_object_desc, 5);
    int flags_int = PyLong_AsLong(flags);
    assert(!ERROR_OCCURRED());

    PyCodeObject *code_object = MAKE_CODEOBJECT(
        filename,
        line_int,
        flags_int,
        function_name,
        argnames,
        NULL, // freevars
        arg_count_int,
        0, // TODO: Missing kw_only_count
        0 // TODO: Missing pos_only_count
    );

    struct Nuitka_FunctionObject *result = Nuitka_Function_New(
        functable___main__[offset],
        code_object->co_name,
#if PYTHON_VERSION >= 0x300
        NULL, // TODO: Not transferring qualname yet
#endif
        code_object,
        defaults,
#if PYTHON_VERSION >= 0x300
        NULL, // kwdefaults are done on the outside currently
        NULL, // TODO: Not transferring annotations
#endif
        module___main__,
        doc,
        NULL,
        0
    );

    return (PyObject *)result;
}

static PyMethodDef _method_def_create_compiled_function = {
    "create_compiled_function",
    (PyCFunction)_create_compiled_function,
    METH_VARARGS | METH_KEYWORDS, NULL
};


#endif

// Internal entry point for module code.
PyObject *modulecode___main__(PyObject *module, struct Nuitka_MetaPathBasedLoaderEntry const *loader_entry) {
    // Report entry to PGO.
    PGO_onModuleEntered("__main__");

    // Store the module for future use.
    module___main__ = module;

    // Modules can be loaded again in case of errors, avoid the init being done again.
    static bool init_done = false;

    if (init_done == false) {
#if defined(_NUITKA_MODULE) && 1
        // In case of an extension module loaded into a process, we need to call
        // initialization here because that's the first and potentially only time
        // we are going called.

        // Initialize the constant values used.
        _initBuiltinModule();
        createGlobalConstants();

        /* Initialize the compiled types of Nuitka. */
        _initCompiledCellType();
        _initCompiledGeneratorType();
        _initCompiledFunctionType();
        _initCompiledMethodType();
        _initCompiledFrameType();

        _initSlotCompare();
#if PYTHON_VERSION >= 0x270
        _initSlotIternext();
#endif

        patchTypeComparison();

        // Enable meta path based loader if not already done.
#ifdef _NUITKA_TRACE
        PRINT_STRING("__main__: Calling setupMetaPathBasedLoader().\n");
#endif
        setupMetaPathBasedLoader();

#if PYTHON_VERSION >= 0x300
        patchInspectModule();
#endif

#endif

        /* The constants only used by this module are created now. */
#ifdef _NUITKA_TRACE
        PRINT_STRING("__main__: Calling createModuleConstants().\n");
#endif
        createModuleConstants();

        /* The code objects used by this module are created now. */
#ifdef _NUITKA_TRACE
        PRINT_STRING("__main__: Calling createModuleCodeObjects().\n");
#endif
        createModuleCodeObjects();

        init_done = true;
    }

    // PRINT_STRING("in init__main__\n");

    moduledict___main__ = MODULE_DICT(module___main__);

#ifdef _NUITKA_PLUGIN_DILL_ENABLED
    registerDillPluginTables(loader_entry->name, &_method_def_reduce_compiled_function, &_method_def_create_compiled_function);
#endif

    // Set "__compiled__" to what version information we have.
    UPDATE_STRING_DICT0(
        moduledict___main__,
        (Nuitka_StringObject *)const_str_plain___compiled__,
        Nuitka_dunder_compiled_value
    );

    // Update "__package__" value to what it ought to be.
    {
#if 1
        UPDATE_STRING_DICT0(
            moduledict___main__,
            (Nuitka_StringObject *)const_str_plain___package__,
            Py_None
        );
#elif 0
        PyObject *module_name = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)const_str_plain___name__);

        UPDATE_STRING_DICT0(
            moduledict___main__,
            (Nuitka_StringObject *)const_str_plain___package__,
            module_name
        );
#else

#if PYTHON_VERSION < 0x300
        PyObject *module_name = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)const_str_plain___name__);
        char const *module_name_cstr = PyString_AS_STRING(module_name);

        char const *last_dot = strrchr(module_name_cstr, '.');

        if (last_dot != NULL) {
            UPDATE_STRING_DICT1(
                moduledict___main__,
                (Nuitka_StringObject *)const_str_plain___package__,
                PyString_FromStringAndSize(module_name_cstr, last_dot - module_name_cstr)
            );
        }
#else
        PyObject *module_name = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)const_str_plain___name__);
        Py_ssize_t dot_index = PyUnicode_Find(module_name, const_str_dot, 0, PyUnicode_GetLength(module_name), -1);

        if (dot_index != -1) {
            UPDATE_STRING_DICT1(
                moduledict___main__,
                (Nuitka_StringObject *)const_str_plain___package__,
                PyUnicode_Substring(module_name, 0, dot_index)
            );
        }
#endif
#endif
    }

    CHECK_OBJECT(module___main__);

    // For deep importing of a module we need to have "__builtins__", so we set
    // it ourselves in the same way than CPython does. Note: This must be done
    // before the frame object is allocated, or else it may fail.

    if (GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)const_str_plain___builtins__) == NULL) {
        PyObject *value = (PyObject *)builtin_module;

        // Check if main module, not a dict then but the module itself.
#if defined(_NUITKA_MODULE) || !1
        value = PyModule_GetDict(value);
#endif

        UPDATE_STRING_DICT0(moduledict___main__, (Nuitka_StringObject *)const_str_plain___builtins__, value);
    }

#if PYTHON_VERSION >= 0x300
    UPDATE_STRING_DICT0(moduledict___main__, (Nuitka_StringObject *)const_str_plain___loader__, (PyObject *)&Nuitka_Loader_Type);
#endif

#if PYTHON_VERSION >= 0x340
// Set the "__spec__" value

#if 1
    // Main modules just get "None" as spec.
    UPDATE_STRING_DICT0(moduledict___main__, (Nuitka_StringObject *)const_str_plain___spec__, Py_None);
#else
    // Other modules get a "ModuleSpec" from the standard mechanism.
    {
        PyObject *bootstrap_module = getImportLibBootstrapModule();
        CHECK_OBJECT(bootstrap_module);

        PyObject *_spec_from_module = PyObject_GetAttrString(bootstrap_module, "_spec_from_module");
        CHECK_OBJECT(_spec_from_module);

        PyObject *spec_value = CALL_FUNCTION_WITH_SINGLE_ARG(_spec_from_module, module___main__);
        Py_DECREF(_spec_from_module);

        // We can assume this to never fail, or else we are in trouble anyway.
        // CHECK_OBJECT(spec_value);

        if (spec_value == NULL) {
            PyErr_PrintEx(0);
            abort();
        }

// Mark the execution in the "__spec__" value.
        SET_ATTRIBUTE(spec_value, const_str_plain__initializing, Py_True);

        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)const_str_plain___spec__, spec_value);
    }
#endif
#endif

    // Temp variables if any
    PyObject *tmp_import_from_1__module = NULL;
    struct Nuitka_FrameObject *frame_a4a3d272fc8190d3932415e3dda5345c;
    NUITKA_MAY_BE_UNUSED char const *type_description_1 = NULL;
    NUITKA_MAY_BE_UNUSED nuitka_void tmp_unused;
    PyObject *exception_type = NULL;
    PyObject *exception_value = NULL;
    PyTracebackObject *exception_tb = NULL;
    NUITKA_MAY_BE_UNUSED int exception_lineno = 0;
    PyObject *exception_keeper_type_1;
    PyObject *exception_keeper_value_1;
    PyTracebackObject *exception_keeper_tb_1;
    NUITKA_MAY_BE_UNUSED int exception_keeper_lineno_1;

    // Module code.
    // Frame without reuse.
    frame_a4a3d272fc8190d3932415e3dda5345c = MAKE_MODULE_FRAME(codeobj_a4a3d272fc8190d3932415e3dda5345c, module___main__);

    // Push the new frame as the currently active one, and we should be exclusively
    // owning it.
    pushFrameStack(frame_a4a3d272fc8190d3932415e3dda5345c);
    assert(Py_REFCNT(frame_a4a3d272fc8190d3932415e3dda5345c) == 2);

    // Framed code:
    {
        PyObject *tmp_imported_value_1;
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 1;
        tmp_imported_value_1 = IMPORT_HARD_SITE();
        if (tmp_imported_value_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 1;

            goto frame_exception_exit_1;
        }
    }
    {
        PyObject *tmp_assign_source_1;
        tmp_assign_source_1 = mod_consts[86];
        UPDATE_STRING_DICT0(moduledict___main__, (Nuitka_StringObject *)mod_consts[87], tmp_assign_source_1);
    }
    {
        PyObject *tmp_assign_source_2;
        tmp_assign_source_2 = mod_consts[88];
        UPDATE_STRING_DICT0(moduledict___main__, (Nuitka_StringObject *)mod_consts[89], tmp_assign_source_2);
    }
    {
        PyObject *tmp_assign_source_3;
        tmp_assign_source_3 = Py_None;
        UPDATE_STRING_DICT0(moduledict___main__, (Nuitka_StringObject *)mod_consts[90], tmp_assign_source_3);
    }
    {
        PyObject *tmp_assign_source_4;
        tmp_assign_source_4 = PyDict_New();
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[91], tmp_assign_source_4);
    }
    {
        PyObject *tmp_assign_source_5;
        PyObject *tmp_name_value_1;
        PyObject *tmp_globals_arg_value_1;
        PyObject *tmp_locals_arg_value_1;
        PyObject *tmp_fromlist_value_1;
        PyObject *tmp_level_value_1;
        tmp_name_value_1 = mod_consts[92];
        tmp_globals_arg_value_1 = (PyObject *)moduledict___main__;
        tmp_locals_arg_value_1 = Py_None;
        tmp_fromlist_value_1 = Py_None;
        tmp_level_value_1 = mod_consts[67];
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 10;
        tmp_assign_source_5 = IMPORT_MODULE5(tmp_name_value_1, tmp_globals_arg_value_1, tmp_locals_arg_value_1, tmp_fromlist_value_1, tmp_level_value_1);
        if (tmp_assign_source_5 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 10;

            goto frame_exception_exit_1;
        }
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[33], tmp_assign_source_5);
    }
    {
        PyObject *tmp_assign_source_6;
        PyObject *tmp_name_value_2;
        PyObject *tmp_globals_arg_value_2;
        PyObject *tmp_locals_arg_value_2;
        PyObject *tmp_fromlist_value_2;
        PyObject *tmp_level_value_2;
        tmp_name_value_2 = mod_consts[93];
        tmp_globals_arg_value_2 = (PyObject *)moduledict___main__;
        tmp_locals_arg_value_2 = Py_None;
        tmp_fromlist_value_2 = Py_None;
        tmp_level_value_2 = mod_consts[67];
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 11;
        tmp_assign_source_6 = IMPORT_MODULE5(tmp_name_value_2, tmp_globals_arg_value_2, tmp_locals_arg_value_2, tmp_fromlist_value_2, tmp_level_value_2);
        if (tmp_assign_source_6 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 11;

            goto frame_exception_exit_1;
        }
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[93], tmp_assign_source_6);
    }
    {
        PyObject *tmp_assign_source_7;
        PyObject *tmp_name_value_3;
        PyObject *tmp_globals_arg_value_3;
        PyObject *tmp_locals_arg_value_3;
        PyObject *tmp_fromlist_value_3;
        PyObject *tmp_level_value_3;
        tmp_name_value_3 = mod_consts[46];
        tmp_globals_arg_value_3 = (PyObject *)moduledict___main__;
        tmp_locals_arg_value_3 = Py_None;
        tmp_fromlist_value_3 = Py_None;
        tmp_level_value_3 = mod_consts[67];
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 12;
        tmp_assign_source_7 = IMPORT_MODULE5(tmp_name_value_3, tmp_globals_arg_value_3, tmp_locals_arg_value_3, tmp_fromlist_value_3, tmp_level_value_3);
        if (tmp_assign_source_7 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 12;

            goto frame_exception_exit_1;
        }
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[46], tmp_assign_source_7);
    }
    {
        PyObject *tmp_assign_source_8;
        PyObject *tmp_name_value_4;
        PyObject *tmp_globals_arg_value_4;
        PyObject *tmp_locals_arg_value_4;
        PyObject *tmp_fromlist_value_4;
        PyObject *tmp_level_value_4;
        tmp_name_value_4 = mod_consts[52];
        tmp_globals_arg_value_4 = (PyObject *)moduledict___main__;
        tmp_locals_arg_value_4 = Py_None;
        tmp_fromlist_value_4 = Py_None;
        tmp_level_value_4 = mod_consts[67];
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 13;
        tmp_assign_source_8 = IMPORT_MODULE5(tmp_name_value_4, tmp_globals_arg_value_4, tmp_locals_arg_value_4, tmp_fromlist_value_4, tmp_level_value_4);
        assert(!(tmp_assign_source_8 == NULL));
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[52], tmp_assign_source_8);
    }
    {
        PyObject *tmp_assign_source_9;
        tmp_assign_source_9 = IMPORT_HARD_OS();
        assert(!(tmp_assign_source_9 == NULL));
        UPDATE_STRING_DICT0(moduledict___main__, (Nuitka_StringObject *)mod_consts[27], tmp_assign_source_9);
    }
    {
        PyObject *tmp_assign_source_10;
        PyObject *tmp_import_name_from_1;
        PyObject *tmp_name_value_5;
        PyObject *tmp_globals_arg_value_5;
        PyObject *tmp_locals_arg_value_5;
        PyObject *tmp_fromlist_value_5;
        PyObject *tmp_level_value_5;
        tmp_name_value_5 = mod_consts[94];
        tmp_globals_arg_value_5 = (PyObject *)moduledict___main__;
        tmp_locals_arg_value_5 = Py_None;
        tmp_fromlist_value_5 = mod_consts[95];
        tmp_level_value_5 = mod_consts[67];
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 15;
        tmp_import_name_from_1 = IMPORT_MODULE5(tmp_name_value_5, tmp_globals_arg_value_5, tmp_locals_arg_value_5, tmp_fromlist_value_5, tmp_level_value_5);
        if (tmp_import_name_from_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 15;

            goto frame_exception_exit_1;
        }
        if (PyModule_Check(tmp_import_name_from_1)) {
            tmp_assign_source_10 = IMPORT_NAME_OR_MODULE(
                tmp_import_name_from_1,
                (PyObject *)moduledict___main__,
                mod_consts[24],
                mod_consts[67]
            );
        } else {
            tmp_assign_source_10 = IMPORT_NAME(tmp_import_name_from_1, mod_consts[24]);
        }

        Py_DECREF(tmp_import_name_from_1);
        if (tmp_assign_source_10 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 15;

            goto frame_exception_exit_1;
        }
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[24], tmp_assign_source_10);
    }
    {
        PyObject *tmp_assign_source_11;
        PyObject *tmp_import_name_from_2;
        PyObject *tmp_name_value_6;
        PyObject *tmp_globals_arg_value_6;
        PyObject *tmp_locals_arg_value_6;
        PyObject *tmp_fromlist_value_6;
        PyObject *tmp_level_value_6;
        tmp_name_value_6 = mod_consts[96];
        tmp_globals_arg_value_6 = (PyObject *)moduledict___main__;
        tmp_locals_arg_value_6 = Py_None;
        tmp_fromlist_value_6 = mod_consts[97];
        tmp_level_value_6 = mod_consts[67];
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 16;
        tmp_import_name_from_2 = IMPORT_MODULE5(tmp_name_value_6, tmp_globals_arg_value_6, tmp_locals_arg_value_6, tmp_fromlist_value_6, tmp_level_value_6);
        if (tmp_import_name_from_2 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 16;

            goto frame_exception_exit_1;
        }
        if (PyModule_Check(tmp_import_name_from_2)) {
            tmp_assign_source_11 = IMPORT_NAME_OR_MODULE(
                tmp_import_name_from_2,
                (PyObject *)moduledict___main__,
                mod_consts[22],
                mod_consts[67]
            );
        } else {
            tmp_assign_source_11 = IMPORT_NAME(tmp_import_name_from_2, mod_consts[22]);
        }

        Py_DECREF(tmp_import_name_from_2);
        if (tmp_assign_source_11 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 16;

            goto frame_exception_exit_1;
        }
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[22], tmp_assign_source_11);
    }
    {
        PyObject *tmp_assign_source_12;
        PyObject *tmp_name_value_7;
        PyObject *tmp_globals_arg_value_7;
        PyObject *tmp_locals_arg_value_7;
        PyObject *tmp_fromlist_value_7;
        PyObject *tmp_level_value_7;
        tmp_name_value_7 = mod_consts[98];
        tmp_globals_arg_value_7 = (PyObject *)moduledict___main__;
        tmp_locals_arg_value_7 = Py_None;
        tmp_fromlist_value_7 = mod_consts[99];
        tmp_level_value_7 = mod_consts[67];
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 18;
        tmp_assign_source_12 = IMPORT_MODULE5(tmp_name_value_7, tmp_globals_arg_value_7, tmp_locals_arg_value_7, tmp_fromlist_value_7, tmp_level_value_7);
        if (tmp_assign_source_12 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 18;

            goto frame_exception_exit_1;
        }
        assert(tmp_import_from_1__module == NULL);
        tmp_import_from_1__module = tmp_assign_source_12;
    }
    // Tried code:
    {
        PyObject *tmp_assign_source_13;
        PyObject *tmp_import_name_from_3;
        CHECK_OBJECT(tmp_import_from_1__module);
        tmp_import_name_from_3 = tmp_import_from_1__module;
        if (PyModule_Check(tmp_import_name_from_3)) {
            tmp_assign_source_13 = IMPORT_NAME_OR_MODULE(
                tmp_import_name_from_3,
                (PyObject *)moduledict___main__,
                mod_consts[0],
                mod_consts[67]
            );
        } else {
            tmp_assign_source_13 = IMPORT_NAME(tmp_import_name_from_3, mod_consts[0]);
        }

        if (tmp_assign_source_13 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 18;

            goto try_except_handler_1;
        }
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[0], tmp_assign_source_13);
    }
    {
        PyObject *tmp_assign_source_14;
        PyObject *tmp_import_name_from_4;
        CHECK_OBJECT(tmp_import_from_1__module);
        tmp_import_name_from_4 = tmp_import_from_1__module;
        if (PyModule_Check(tmp_import_name_from_4)) {
            tmp_assign_source_14 = IMPORT_NAME_OR_MODULE(
                tmp_import_name_from_4,
                (PyObject *)moduledict___main__,
                mod_consts[1],
                mod_consts[67]
            );
        } else {
            tmp_assign_source_14 = IMPORT_NAME(tmp_import_name_from_4, mod_consts[1]);
        }

        if (tmp_assign_source_14 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 18;

            goto try_except_handler_1;
        }
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[1], tmp_assign_source_14);
    }
    goto try_end_1;
    // Exception handler code:
    try_except_handler_1:;
    exception_keeper_type_1 = exception_type;
    exception_keeper_value_1 = exception_value;
    exception_keeper_tb_1 = exception_tb;
    exception_keeper_lineno_1 = exception_lineno;
    exception_type = NULL;
    exception_value = NULL;
    exception_tb = NULL;
    exception_lineno = 0;

    CHECK_OBJECT(tmp_import_from_1__module);
    Py_DECREF(tmp_import_from_1__module);
    tmp_import_from_1__module = NULL;
    // Re-raise.
    exception_type = exception_keeper_type_1;
    exception_value = exception_keeper_value_1;
    exception_tb = exception_keeper_tb_1;
    exception_lineno = exception_keeper_lineno_1;

    goto frame_exception_exit_1;
    // End of try:
    try_end_1:;
    CHECK_OBJECT(tmp_import_from_1__module);
    Py_DECREF(tmp_import_from_1__module);
    tmp_import_from_1__module = NULL;
    {
        PyObject *tmp_called_instance_1;
        PyObject *tmp_call_result_1;
        tmp_called_instance_1 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[93]);

        if (unlikely(tmp_called_instance_1 == NULL)) {
            tmp_called_instance_1 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[93]);
        }

        if (tmp_called_instance_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 20;

            goto frame_exception_exit_1;
        }
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 20;
        tmp_call_result_1 = CALL_METHOD_WITH_SINGLE_ARG(
            tmp_called_instance_1,
            mod_consts[100],
            PyTuple_GET_ITEM(mod_consts[101], 0)
        );

        if (tmp_call_result_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 20;

            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_1);
    }
    {
        PyObject *tmp_assign_source_15;
        PyObject *tmp_called_value_1;
        PyObject *tmp_expression_value_1;
        PyObject *tmp_expression_value_2;
        PyObject *tmp_args_element_value_1;
        PyObject *tmp_called_value_2;
        PyObject *tmp_expression_value_3;
        PyObject *tmp_expression_value_4;
        PyObject *tmp_args_element_value_2;
        tmp_expression_value_2 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[27]);

        if (unlikely(tmp_expression_value_2 == NULL)) {
            tmp_expression_value_2 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[27]);
        }

        if (tmp_expression_value_2 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 21;

            goto frame_exception_exit_1;
        }
        tmp_expression_value_1 = LOOKUP_ATTRIBUTE(tmp_expression_value_2, mod_consts[28]);
        if (tmp_expression_value_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 21;

            goto frame_exception_exit_1;
        }
        tmp_called_value_1 = LOOKUP_ATTRIBUTE(tmp_expression_value_1, mod_consts[102]);
        Py_DECREF(tmp_expression_value_1);
        if (tmp_called_value_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 21;

            goto frame_exception_exit_1;
        }
        tmp_expression_value_4 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[27]);

        if (unlikely(tmp_expression_value_4 == NULL)) {
            tmp_expression_value_4 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[27]);
        }

        if (tmp_expression_value_4 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_1);

            exception_lineno = 21;

            goto frame_exception_exit_1;
        }
        tmp_expression_value_3 = LOOKUP_ATTRIBUTE(tmp_expression_value_4, mod_consts[28]);
        if (tmp_expression_value_3 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_1);

            exception_lineno = 21;

            goto frame_exception_exit_1;
        }
        tmp_called_value_2 = LOOKUP_ATTRIBUTE(tmp_expression_value_3, mod_consts[103]);
        Py_DECREF(tmp_expression_value_3);
        if (tmp_called_value_2 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_1);

            exception_lineno = 21;

            goto frame_exception_exit_1;
        }
        tmp_args_element_value_2 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[89]);

        if (unlikely(tmp_args_element_value_2 == NULL)) {
            tmp_args_element_value_2 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[89]);
        }

        if (tmp_args_element_value_2 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_1);
            Py_DECREF(tmp_called_value_2);

            exception_lineno = 21;

            goto frame_exception_exit_1;
        }
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 21;
        tmp_args_element_value_1 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_2, tmp_args_element_value_2);
        Py_DECREF(tmp_called_value_2);
        if (tmp_args_element_value_1 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);
            Py_DECREF(tmp_called_value_1);

            exception_lineno = 21;

            goto frame_exception_exit_1;
        }
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 21;
        tmp_assign_source_15 = CALL_FUNCTION_WITH_SINGLE_ARG(tmp_called_value_1, tmp_args_element_value_1);
        Py_DECREF(tmp_called_value_1);
        Py_DECREF(tmp_args_element_value_1);
        if (tmp_assign_source_15 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 21;

            goto frame_exception_exit_1;
        }
        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[23], tmp_assign_source_15);
    }
    {
        PyObject *tmp_assign_source_16;


        tmp_assign_source_16 = MAKE_FUNCTION___main__$$$function__1_main();

        UPDATE_STRING_DICT1(moduledict___main__, (Nuitka_StringObject *)mod_consts[104], tmp_assign_source_16);
    }
    {
        PyObject *tmp_called_value_3;
        PyObject *tmp_call_result_2;
        tmp_called_value_3 = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)mod_consts[104]);

        if (unlikely(tmp_called_value_3 == NULL)) {
            tmp_called_value_3 = GET_MODULE_VARIABLE_VALUE_FALLBACK(mod_consts[104]);
        }

        assert(!(tmp_called_value_3 == NULL));
        frame_a4a3d272fc8190d3932415e3dda5345c->m_frame.f_lineno = 97;
        tmp_call_result_2 = CALL_FUNCTION_NO_ARGS(tmp_called_value_3);
        if (tmp_call_result_2 == NULL) {
            assert(ERROR_OCCURRED());

            FETCH_ERROR_OCCURRED(&exception_type, &exception_value, &exception_tb);


            exception_lineno = 97;

            goto frame_exception_exit_1;
        }
        Py_DECREF(tmp_call_result_2);
    }

    // Restore frame exception if necessary.
#if 0
    RESTORE_FRAME_EXCEPTION(frame_a4a3d272fc8190d3932415e3dda5345c);
#endif
    popFrameStack();

    assertFrameObject(frame_a4a3d272fc8190d3932415e3dda5345c);

    goto frame_no_exception_1;

    frame_exception_exit_1:;
#if 0
    RESTORE_FRAME_EXCEPTION(frame_a4a3d272fc8190d3932415e3dda5345c);
#endif

    if (exception_tb == NULL) {
        exception_tb = MAKE_TRACEBACK(frame_a4a3d272fc8190d3932415e3dda5345c, exception_lineno);
    } else if (exception_tb->tb_frame != &frame_a4a3d272fc8190d3932415e3dda5345c->m_frame) {
        exception_tb = ADD_TRACEBACK(exception_tb, frame_a4a3d272fc8190d3932415e3dda5345c, exception_lineno);
    }

    // Put the previous frame back on top.
    popFrameStack();

    // Return the error.
    goto module_exception_exit;

    frame_no_exception_1:;

    // Report to PGO about leaving the module without error.
    PGO_onModuleExit("__main__", false);

    return module___main__;
    module_exception_exit:

#if defined(_NUITKA_MODULE) && 1
    {
        PyObject *module_name = GET_STRING_DICT_VALUE(moduledict___main__, (Nuitka_StringObject *)const_str_plain___name__);

        if (module_name != NULL) {
            Nuitka_DelModule(module_name);
        }
    }
#endif
    PGO_onModuleExit("__main__", false);

    RESTORE_ERROR_OCCURRED(exception_type, exception_value, exception_tb);
    return NULL;
}
