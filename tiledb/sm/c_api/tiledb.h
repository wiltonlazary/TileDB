/**
 * @file   tiledb.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2021 TileDB, Inc.
 * @copyright Copyright (c) 2016 MIT and Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * @section DESCRIPTION
 *
 * This file declares the C API for TileDB.
 */

#ifndef TILEDB_H
#define TILEDB_H

#include "tiledb_version.h"

#ifdef TILEDB_NO_API_DEPRECATION_WARNINGS
// Define these before including tiledb_export.h to avoid their normal
// definitions.
#ifndef TILEDB_DEPRECATED
#define TILEDB_DEPRECATED
#endif
#ifndef TILEDB_DEPRECATED_EXPORT
#define TILEDB_DEPRECATED_EXPORT
#endif
#endif

/*
 * Common definitions for export, noexcept, etc.
 */
#include "tiledb/api/c_api/api_external_common.h"

/*
 * API sections
 */
#include "tiledb/api/c_api/config/config_api_external.h"
#include "tiledb/api/c_api/context/context_api_external.h"
#include "tiledb/api/c_api/error/error_api_external.h"
#include "tiledb/api/c_api/filesystem/filesystem_api_external.h"
#include "tiledb/api/c_api/filter/filter_api_external.h"
#include "tiledb/api/c_api/filter_list/filter_list_api_external.h"

#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ****************************** */
/*          TILEDB ENUMS          */
/* ****************************** */

/** TileDB object type. */
typedef enum {
/** Helper macro for defining object type enums. */
#define TILEDB_OBJECT_TYPE_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_OBJECT_TYPE_ENUM
} tiledb_object_t;

/** TileDB query type. */
typedef enum {
/** Helper macro for defining query type enums. */
#define TILEDB_QUERY_TYPE_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_QUERY_TYPE_ENUM
} tiledb_query_type_t;

/** Query status. */
typedef enum {
/** Helper macro for defining query status enums. */
#define TILEDB_QUERY_STATUS_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_QUERY_STATUS_ENUM
} tiledb_query_status_t;

/** Query condition operator. */
typedef enum {
/** Helper macro for defining query condition operator enums. */
#define TILEDB_QUERY_CONDITION_OP_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_QUERY_CONDITION_OP_ENUM
} tiledb_query_condition_op_t;

/** Query condition combination operator. */
typedef enum {
/** Helper macro for defining query condition combination operator enums. */
#define TILEDB_QUERY_CONDITION_COMBINATION_OP_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_QUERY_CONDITION_COMBINATION_OP_ENUM
} tiledb_query_condition_combination_op_t;

/** TileDB datatype. */
typedef enum {
/** Helper macro for defining datatype enums. */
#define TILEDB_DATATYPE_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_DATATYPE_ENUM
#ifdef TILEDB_CHAR
#def TILEDB_CHAR_VAL TILEDB_CHAR
#undef TILEDB_CHAR
#define TILEDB_CHAR TILEDB_DEPRECATED TILEDB_CHAR_VAL
#undef TILEDB_CHAR_VAL
#endif
#ifdef TILEDB_STRING_UCS2
#def TILEDB_STRING_UCS2_VAL TILEDB_STRING_UCS2
#undef TILEDB_STRING_UCS2
#define TILEDB_STRING_UCS2 TILEDB_DEPRECATED TILEDB_STRING_UCS2_VAL
#undef TILEDB_STRING_UCS2_VAL
#endif
#ifdef TILEDB_STRING_UCS4
#def TILEDB_STRING_UCS4_VAL TILEDB_STRING_UCS4
#undef TILEDB_STRING_UCS4
#define TILEDB_STRING_UCS4 TILEDB_DEPRECATED TILEDB_STRING_UCS4_VAL
#undef TILEDB_STRING_UCS4_VAL
#endif
#ifdef TILEDB_ANY
#def TILEDB_ANY_VAL TILEDB_ANY
#undef TILEDB_ANY
#define TILEDB_ANY TILEDB_DEPRECATED TILEDB_ANY_VAL
#undef TILEDB_ANY_VAL
#endif
} tiledb_datatype_t;

/** Array type. */
typedef enum {
/** Helper macro for defining array type enums. */
#define TILEDB_ARRAY_TYPE_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_ARRAY_TYPE_ENUM
} tiledb_array_type_t;

/** Tile or cell layout. */
typedef enum {
/** Helper macro for defining layout type enums. */
#define TILEDB_LAYOUT_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_LAYOUT_ENUM
} tiledb_layout_t;

/** Encryption type. */
typedef enum {
/** Helper macro for defining encryption enums. */
#define TILEDB_ENCRYPTION_TYPE_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_ENCRYPTION_TYPE_ENUM
} tiledb_encryption_type_t;

/** Walk traversal order. */
typedef enum {
/** Helper macro for defining walk order enums. */
#define TILEDB_WALK_ORDER_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_WALK_ORDER_ENUM
} tiledb_walk_order_t;

/** VFS mode. */
typedef enum {
/** Helper macro for defining VFS mode enums. */
#define TILEDB_VFS_MODE_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_VFS_MODE_ENUM
} tiledb_vfs_mode_t;

/** MIME Type. */
typedef enum {
/** Helper macro for defining MimeType enums. */
#define TILEDB_MIME_TYPE_ENUM(id) TILEDB_##id
#include "tiledb_enum.h"
#undef TILEDB_MIME_TYPE_ENUM
} tiledb_mime_type_t;

/* ****************************** */
/*       ENUMS TO/FROM STR        */
/* ****************************** */

/**
 * Returns a string representation of the given query type.
 *
 * @param query_type Query type
 * @param str Set to point to a constant string representation of the query type
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_type_to_str(
    tiledb_query_type_t query_type, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a query type from the given string.
 *
 * @param str String representation to parse
 * @param query_type Set to the parsed query type
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_type_from_str(
    const char* str, tiledb_query_type_t* query_type) TILEDB_NOEXCEPT;

/**
 * Returns a string representation of the given object type.
 *
 * @param object_type Object type
 * @param str Set to point to a constant string representation of the object
 * type
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_object_type_to_str(
    tiledb_object_t object_type, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a object type from the given string.
 *
 * @param str String representation to parse
 * @param object_type Set to the parsed object type
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_object_type_from_str(
    const char* str, tiledb_object_t* object_type) TILEDB_NOEXCEPT;

/**
 * Returns a string representation of the given filesystem.
 *
 * @param filesystem Filesystem
 * @param str Set to point to a constant string representation of the filesystem
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_filesystem_to_str(
    tiledb_filesystem_t filesystem, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a filesystem from the given string.
 *
 * @param str String representation to parse
 * @param filesystem Set to the parsed filesystem
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_filesystem_from_str(
    const char* str, tiledb_filesystem_t* filesystem) TILEDB_NOEXCEPT;

/**
 * Returns a string representation of the given datatype.
 *
 * @param datatype Datatype
 * @param str Set to point to a constant string representation of the datatype
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_datatype_to_str(
    tiledb_datatype_t datatype, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a datatype from the given string.
 *
 * @param str String representation to parse
 * @param datatype Set to the parsed datatype
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_datatype_from_str(
    const char* str, tiledb_datatype_t* datatype) TILEDB_NOEXCEPT;

/**
 * Returns a string representation of the given array type.
 *
 * @param array_type Array type
 * @param str Set to point to a constant string representation of the array type
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_type_to_str(
    tiledb_array_type_t array_type, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a array type from the given string.
 *
 * @param str String representation to parse
 * @param array_type Set to the parsed array type
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_type_from_str(
    const char* str, tiledb_array_type_t* array_type) TILEDB_NOEXCEPT;

/**
 * Returns a string representation of the given layout.
 *
 * @param layout Layout
 * @param str Set to point to a constant string representation of the layout
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t
tiledb_layout_to_str(tiledb_layout_t layout, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a layout from the given string.
 *
 * @param str String representation to parse
 * @param layout Set to the parsed layout
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_layout_from_str(
    const char* str, tiledb_layout_t* layout) TILEDB_NOEXCEPT;

/**
 * Returns a string representation of the given encryption type.
 *
 * @param encryption_type Encryption type
 * @param str Set to point to a constant string representation of the encryption
 * type
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_encryption_type_to_str(
    tiledb_encryption_type_t encryption_type, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a encryption type from the given string.
 *
 * @param str String representation to parse
 * @param encryption_type Set to the parsed encryption type
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_encryption_type_from_str(
    const char* str, tiledb_encryption_type_t* encryption_type) TILEDB_NOEXCEPT;

/**
 * Returns a string representation of the given query status.
 *
 * @param query_status Query status
 * @param str Set to point to a constant string representation of the query
 * status
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_status_to_str(
    tiledb_query_status_t query_status, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a query status from the given string.
 *
 * @param str String representation to parse
 * @param query_status Set to the parsed query status
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_status_from_str(
    const char* str, tiledb_query_status_t* query_status) TILEDB_NOEXCEPT;

/**
 * Returns a string representation of the given walk order.
 *
 * @param walk_order Walk order
 * @param str Set to point to a constant string representation of the walk order
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_walk_order_to_str(
    tiledb_walk_order_t walk_order, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a walk order from the given string.
 *
 * @param str String representation to parse
 * @param walk_order Set to the parsed walk order
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_walk_order_from_str(
    const char* str, tiledb_walk_order_t* walk_order) TILEDB_NOEXCEPT;

/**
 * Returns a string representation of the given VFS mode.
 *
 * @param vfs_mode VFS mode
 * @param str Set to point to a constant string representation of the VFS mode
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_mode_to_str(
    tiledb_vfs_mode_t vfs_mode, const char** str) TILEDB_NOEXCEPT;

/**
 * Parses a VFS mode from the given string.
 *
 * @param str String representation to parse
 * @param vfs_mode Set to the parsed VFS mode
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_mode_from_str(
    const char* str, tiledb_vfs_mode_t* vfs_mode) TILEDB_NOEXCEPT;

/* ****************************** */
/*            CONSTANTS           */
/* ****************************** */

/**
 * Returns a special name indicating the coordinates attribute.
 *
 * The coordinate buffer has been deprecated. Set the coordinates for
 * each individual dimension with the `set_buffer` API. Consult the current
 * documentation for more information.
 */
TILEDB_DEPRECATED_EXPORT const char* tiledb_coords(void) TILEDB_NOEXCEPT;

/** Returns a special value indicating a variable number of elements. */
TILEDB_EXPORT uint32_t tiledb_var_num(void) TILEDB_NOEXCEPT;

/** Returns the maximum path length on the current platform. */
TILEDB_EXPORT uint32_t tiledb_max_path(void) TILEDB_NOEXCEPT;

/**
 * Returns the size (in bytes) of an offset (used in variable-sized
 * attributes).
 */
TILEDB_EXPORT uint64_t tiledb_offset_size(void) TILEDB_NOEXCEPT;

/**
 * Returns the input datatype size for a given type. Returns zero if the type is
 * not valid.
 */
TILEDB_EXPORT uint64_t tiledb_datatype_size(tiledb_datatype_t type)
    TILEDB_NOEXCEPT;

/** Returns the current time in milliseconds. */
TILEDB_EXPORT uint64_t tiledb_timestamp_now_ms(void) TILEDB_NOEXCEPT;

/** Returns a special name indicating the timestamps attribute. */
TILEDB_EXPORT const char* tiledb_timestamps(void) TILEDB_NOEXCEPT;

/**
 * @name Constants wrapping special functions
 */
/**@{*/
/** A special name indicating the coordinates attribute. */
#define TILEDB_COORDS tiledb_coords()
/** A special value indicating a variable number of elements. */
#define TILEDB_VAR_NUM tiledb_var_num()
/** The maximum path length on the current platform. */
#define TILEDB_MAX_PATH tiledb_max_path()
/** The size (in bytes) of an offset (used in variable-sized attributes). */
#define TILEDB_OFFSET_SIZE tiledb_offset_size()
/** The current time in milliseconds. */
#define TILEDB_TIMESTAMP_NOW_MS tiledb_timestamp_now_ms()
/** A special name indicating the timestamps attribute. */
#define TILEDB_TIMESTAMPS tiledb_timestamps()
/**@}*/

/* ****************************** */
/*            VERSION             */
/* ****************************** */

/**
 *  Retrieves the version of the TileDB library currently being used.
 *
 *  @param major Will store the major version number.
 *  @param minor Will store the minor version number.
 *  @param rev Will store the revision (patch) number.
 */
TILEDB_EXPORT void tiledb_version(int32_t* major, int32_t* minor, int32_t* rev)
    TILEDB_NOEXCEPT;

/* ********************************* */
/*           TILEDB TYPES            */
/* ********************************* */

/** An array object. */
typedef struct tiledb_array_t tiledb_array_t;

/** A subarray object. */
typedef struct tiledb_subarray_t tiledb_subarray_t;

/** A generic buffer object. */
typedef struct tiledb_buffer_t tiledb_buffer_t;

/** A generic buffer list object. */
typedef struct tiledb_buffer_list_t tiledb_buffer_list_t;

/** A TileDB attribute. */
typedef struct tiledb_attribute_t tiledb_attribute_t;

/** A TileDB array schema. */
typedef struct tiledb_array_schema_t tiledb_array_schema_t;

/** A TileDB dimension. */
typedef struct tiledb_dimension_t tiledb_dimension_t;

/** A TileDB domain. */
typedef struct tiledb_domain_t tiledb_domain_t;

/** A TileDB query. */
typedef struct tiledb_query_t tiledb_query_t;

/** A TileDB query condition object. */
typedef struct tiledb_query_condition_t tiledb_query_condition_t;

/** A virtual filesystem object. */
typedef struct tiledb_vfs_t tiledb_vfs_t;

/** A virtual filesystem file handle. */
typedef struct tiledb_vfs_fh_t tiledb_vfs_fh_t;

/** A fragment info object. */
typedef struct tiledb_fragment_info_t tiledb_fragment_info_t;

/** A group object. */
typedef struct tiledb_group_t tiledb_group_t;

/* ********************************* */
/*              BUFFER               */
/* ********************************* */

/**
 * Creates an empty buffer object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_t* buffer;
 * tiledb_buffer_alloc(ctx, &buffer);
 * @endcode
 *
 * @param ctx TileDB context
 * @param buffer The buffer to be created
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_alloc(
    tiledb_ctx_t* ctx, tiledb_buffer_t** buffer) TILEDB_NOEXCEPT;

/**
 * Destroys a TileDB buffer, freeing associated memory.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_t* buffer;
 * tiledb_buffer_alloc(ctx, &buffer);
 * tiledb_buffer_free(&buffer);
 * @endcode
 *
 * @param buffer The buffer to be destroyed.
 */
TILEDB_EXPORT void tiledb_buffer_free(tiledb_buffer_t** buffer) TILEDB_NOEXCEPT;

/**
 * Sets a datatype for the given buffer. The default datatype is `TILEDB_UINT8`.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_t* buffer;
 * tiledb_buffer_alloc(ctx, &buffer);
 * tiledb_buffer_set_type(ctx, buffer, TILEDB_INT32);
 * @endcode
 *
 * @param ctx TileDB context
 * @param buffer TileDB buffer instance
 * @param datatype The datatype to set on the buffer.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_set_type(
    tiledb_ctx_t* ctx,
    tiledb_buffer_t* buffer,
    tiledb_datatype_t datatype) TILEDB_NOEXCEPT;

/**
 * Gets the datatype from the given buffer.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_datatype_t type;
 * tiledb_buffer_get_type(ctx, buffer, &type);
 * @endcode
 *
 * @param ctx TileDB context
 * @param buffer TileDB buffer instance
 * @param datatype Set to the datatype of the buffer.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_get_type(
    tiledb_ctx_t* ctx,
    const tiledb_buffer_t* buffer,
    tiledb_datatype_t* datatype) TILEDB_NOEXCEPT;

/**
 * Gets a pointer to the current allocation and the current number of bytes in
 * the specified buffer object.
 *
 * @note For string buffers allocated by TileDB, the number of bytes includes
 * the terminating NULL byte.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_t* buffer;
 * tiledb_buffer_alloc(ctx, &buffer);
 * void* data;
 * uint64_t num_bytes;
 * tiledb_buffer_get_data(ctx, buffer, &data, num_bytes);
 * // data == NULL and num_bytes == 0 because the buffer is currently empty.
 * tiledb_buffer_free(&buffer);
 * @endcode
 *
 * @param ctx TileDB context
 * @param buffer TileDB buffer instance
 * @param data The pointer to the data to be retrieved.
 * @param num_bytes Set to the number of bytes in the buffer.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_get_data(
    tiledb_ctx_t* ctx,
    const tiledb_buffer_t* buffer,
    void** data,
    uint64_t* num_bytes) TILEDB_NOEXCEPT;

/**
 * Sets (wraps) a pre-allocated region of memory with the given buffer object.
 * This does not perform a copy.
 *
 * @note The TileDB buffer object does not take ownership of the allocation
 * set with this function. That means the call to `tiledb_buffer_free` will not
 * free a user allocation set via `tiledb_buffer_set_buffer`.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_t* buffer;
 * tiledb_buffer_alloc(ctx, &buffer);
 *
 * void* my_data = malloc(100);
 * tiledb_buffer_set_data(ctx, buffer, my_data, 100);
 *
 * void* data;
 * uint64_t num_bytes;
 * tiledb_buffer_get_data(ctx, buffer, &data, num_bytes);
 * assert(data == my_data);
 * assert(num_bytes == 100);
 *
 * tiledb_buffer_free(&buffer);
 * free(my_data);
 * @endcode
 *
 * @param ctx TileDB context
 * @param buffer TileDB buffer instance
 * @param data Pre-allocated region of memory to wrap with this buffer.
 * @param size Size (in bytes) of the region pointed to by data.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_set_data(
    tiledb_ctx_t* ctx, tiledb_buffer_t* buffer, void* data, uint64_t size);

/* ********************************* */
/*            BUFFER LIST            */
/* ********************************* */

/**
 * Creates an empty buffer list object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_list_t* buffer_list;
 * tiledb_buffer_list_alloc(ctx, &buffer_list);
 * @endcode
 *
 * @param ctx TileDB context
 * @param buffer_list The buffer list to be created
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_list_alloc(
    tiledb_ctx_t* ctx, tiledb_buffer_list_t** buffer_list) TILEDB_NOEXCEPT;

/**
 * Destroys a TileDB buffer list, freeing associated memory.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_t* buffer_list;
 * tiledb_buffer_list_alloc(ctx, &buffer_list);
 * tiledb_buffer_list_free(&buffer_list);
 * @endcode
 *
 * @param buffer_list The buffer list to be destroyed.
 */
TILEDB_EXPORT void tiledb_buffer_list_free(tiledb_buffer_list_t** buffer_list)
    TILEDB_NOEXCEPT;

/**
 * Gets the number of buffers in the buffer list.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_list_t* buffer_list;
 * tiledb_buffer_list_alloc(ctx, &buffer_list);
 * uint64_t num_buffers;
 * tiledb_buffer_list_get_num_buffers(ctx, buffer_list, &num_buffers);
 * // num_buffers == 0 because the list is empty.
 * @endcode
 *
 * @param ctx TileDB context.
 * @param buffer_list The buffer list.
 * @param num_buffers Set to the number of buffers in the buffer list.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_list_get_num_buffers(
    tiledb_ctx_t* ctx,
    const tiledb_buffer_list_t* buffer_list,
    uint64_t* num_buffers) TILEDB_NOEXCEPT;

/**
 * Gets the buffer at the given index in the buffer list. The returned buffer
 * object is simply a pointer to memory managed by the underlying buffer
 * list, meaning this function does not perform a copy.
 *
 * It is the caller's responsibility to free the returned buffer with
 * `tiledb_buffer_free`. Since the returned buffer object does not "own" the
 * underlying allocation, the underlying allocation is not freed when freeing it
 * with `tiledb_buffer_free`.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_list_t* buffer_list;
 * // Create and populate the buffer_list
 *
 * // Get the buffer at index 0.
 * tiledb_buffer_t *buff0;
 * tiledb_buffer_list_get_buffer(ctx, buffer_list, 0, &buff0);
 *
 * // Always free the returned buffer object
 * tiledb_buffer_free(&buff0);
 * tiledb_buffer_list_free(&buffer_list);
 * @endcode
 *
 * @param ctx TileDB context.
 * @param buffer_list The buffer list.
 * @param buffer_idx Index of buffer to get from the buffer list.
 * @param buffer Set to a newly allocated buffer object pointing to the
 *    underlying allocation in the buffer list corresponding to the buffer.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_list_get_buffer(
    tiledb_ctx_t* ctx,
    const tiledb_buffer_list_t* buffer_list,
    uint64_t buffer_idx,
    tiledb_buffer_t** buffer) TILEDB_NOEXCEPT;

/**
 * Gets the total number of bytes in the buffers in the buffer list.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_list_t* buffer_list;
 * tiledb_buffer_list_alloc(ctx, &buffer_list);
 * uint64_t total_size;
 * tiledb_buffer_list_get_total_size(ctx, buffer_list, &total_size);
 * // total_size == 0 because the list is empty.
 * @endcode
 *
 * @param ctx TileDB context.
 * @param buffer_list The buffer list.
 * @param total_size Set to the total number of bytes in the buffers in the
 *    buffer list.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_list_get_total_size(
    tiledb_ctx_t* ctx,
    const tiledb_buffer_list_t* buffer_list,
    uint64_t* total_size) TILEDB_NOEXCEPT;

/**
 * Copies and concatenates all the data in the buffer list into a new buffer.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_buffer_t* buff;
 * tiledb_buffer_list_flatten(ctx, buffer_list, &buff);
 * // ...
 * tiledb_buffer_free(&buff);
 * @endcode
 *
 * @param ctx TileDB context.
 * @param buffer_list The buffer list.
 * @param buffer Will be set to a newly allocated buffer holding a copy of the
 *    concatenated data from the buffer list.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_buffer_list_flatten(
    tiledb_ctx_t* ctx,
    const tiledb_buffer_list_t* buffer_list,
    tiledb_buffer_t** buffer) TILEDB_NOEXCEPT;

/* ********************************* */
/*                GROUP              */
/* ********************************* */

/**
 * Creates a new TileDB group.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_group_create(ctx, "my_group");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param group_uri The group URI.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t
tiledb_group_create(tiledb_ctx_t* ctx, const char* group_uri) TILEDB_NOEXCEPT;

/* ********************************* */
/*            ATTRIBUTE              */
/* ********************************* */

/**
 * Creates a TileDB attribute.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_attribute_t* attr;
 * tiledb_attribute_alloc(ctx, "my_attr", TILEDB_INT32, &attr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param name The attribute name. Providing an empty string for the name
 * creates an anonymous attribute.
 * @param type The attribute type.
 * @param attr The TileDB attribute to be created.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 *
 * @note The default number of values per cell is 1.
 */
TILEDB_EXPORT int32_t tiledb_attribute_alloc(
    tiledb_ctx_t* ctx,
    const char* name,
    tiledb_datatype_t type,
    tiledb_attribute_t** attr) TILEDB_NOEXCEPT;

/**
 * Destroys a TileDB attribute, freeing associated memory.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_attribute_t* attr;
 * tiledb_attribute_alloc(ctx, "my_attr", TILEDB_INT32, &attr);
 * tiledb_attribute_free(&attr);
 * @endcode
 *
 * @param attr The attribute to be destroyed.
 */
TILEDB_EXPORT void tiledb_attribute_free(tiledb_attribute_t** attr)
    TILEDB_NOEXCEPT;

/**
 * Sets the nullability of an attribute.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_attribute_set_nullable(ctx, attr, 1);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The target attribute.
 * @param nullable Non-zero if the attribute is nullable.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_set_nullable(
    tiledb_ctx_t* ctx,
    tiledb_attribute_t* attr,
    uint8_t nullable) TILEDB_NOEXCEPT;

/**
 * Sets the filter list for an attribute.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_filter_list_alloc(ctx, &filter_list);
 * tiledb_filter_list_add_filter(ctx, filter_list, filter);
 * tiledb_attribute_set_filter_list(ctx, attr, filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The target attribute.
 * @param filter_list The filter_list to be set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_set_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_attribute_t* attr,
    tiledb_filter_list_t* filter_list) TILEDB_NOEXCEPT;

/**
 * Sets the number of values per cell for an attribute. If this is not
 * used, the default is `1`.
 *
 * **Examples:**
 *
 * For a fixed-sized attribute:
 *
 * @code{.c}
 * tiledb_attribute_set_cell_val_num(ctx, attr, 3);
 * @endcode
 *
 * For a variable-sized attribute:
 *
 * @code{.c}
 * tiledb_attribute_set_cell_val_num(ctx, attr, TILEDB_VAR_NUM);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The target attribute.
 * @param cell_val_num The number of values per cell.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_set_cell_val_num(
    tiledb_ctx_t* ctx,
    tiledb_attribute_t* attr,
    uint32_t cell_val_num) TILEDB_NOEXCEPT;

/**
 * Retrieves the attribute name.
 *
 * **Example:**
 *
 * @code{.c}
 * const char* attr_name;
 * tiledb_attribute_get_name(ctx, attr, &attr_name);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The attribute.
 * @param name The name to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_get_name(
    tiledb_ctx_t* ctx,
    const tiledb_attribute_t* attr,
    const char** name) TILEDB_NOEXCEPT;

/**
 * Retrieves the attribute type.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_datatype_t attr_type;
 * tiledb_attribute_get_type(ctx, attr, &attr_type);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The attribute.
 * @param type The type to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_get_type(
    tiledb_ctx_t* ctx,
    const tiledb_attribute_t* attr,
    tiledb_datatype_t* type) TILEDB_NOEXCEPT;

/**
 * Sets the nullability of an attribute.
 *
 * **Example:**
 *
 * @code{.c}
 * uint8_t nullable;
 * tiledb_attribute_get_nullable(ctx, attr, &nullable);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The target attribute.
 * @param nullable Output argument, non-zero for nullable and zero
 *    for non-nullable.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_get_nullable(
    tiledb_ctx_t* ctx,
    tiledb_attribute_t* attr,
    uint8_t* nullable) TILEDB_NOEXCEPT;

/**
 * Retrieves the filter list for an attribute.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_attribute_get_filter_list(ctx, attr, &filter_list);
 * tiledb_filter_list_free(&filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The target attribute.
 * @param filter_list The filter list to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_get_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_attribute_t* attr,
    tiledb_filter_list_t** filter_list) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of values per cell for the attribute. For variable-sized
 * attributes result is TILEDB_VAR_NUM.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t num;
 * tiledb_attribute_get_cell_val_num(ctx, attr, &num);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The attribute.
 * @param cell_val_num The number of values per cell to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_get_cell_val_num(
    tiledb_ctx_t* ctx,
    const tiledb_attribute_t* attr,
    uint32_t* cell_val_num) TILEDB_NOEXCEPT;

/**
 * Retrieves the cell size for this attribute.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t cell_size;
 * tiledb_attribute_get_cell_size(ctx, attr, &cell_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The attribute.
 * @param cell_size The cell size to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_get_cell_size(
    tiledb_ctx_t* ctx,
    const tiledb_attribute_t* attr,
    uint64_t* cell_size) TILEDB_NOEXCEPT;

/**
 * Dumps the contents of an attribute in ASCII form to some output (e.g.,
 * file or stdout).
 *
 * **Example:**
 *
 * The following prints the attribute dump to standard output.
 *
 * @code{.c}
 * tiledb_attribute_dump(ctx, attr, stdout);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The attribute.
 * @param out The output.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error./
 */
TILEDB_EXPORT int32_t tiledb_attribute_dump(
    tiledb_ctx_t* ctx,
    const tiledb_attribute_t* attr,
    FILE* out) TILEDB_NOEXCEPT;

/**
 * Sets the default fill value for the input attribute. This value will
 * be used for the input attribute whenever querying (1) an empty cell in
 * a dense array, or (2) a non-empty cell (in either dense or sparse array)
 * when values on the input attribute are missing (e.g., if the user writes
 * a subset of the attributes in a write operation).
 *
 * Applicable to var-sized attributes.
 *
 * **Example:**
 *
 * @code{.c}
 * // Assumming a int32 attribute
 * int32_t value = 0;
 * uint64_t size = sizeof(value);
 * tiledb_attribute_set_fill_value(ctx, attr, &value, size);
 *
 * // Assumming a var char attribute
 * const char* value = "foo";
 * uint64_t size = strlen(value);
 * tiledb_attribute_set_fill_value(ctx, attr, value, size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The target attribute.
 * @param value The fill value to set.
 * @param size The fill value size in bytes.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note A call to `tiledb_attribute_cell_val_num` sets the fill value
 *     of the attribute to its default. Therefore, make sure you invoke
 *     `tiledb_attribute_set_fill_value` after deciding on the number
 *     of values this attribute will hold in each cell.
 *
 * @note For fixed-sized attributes, the input `size` should be equal
 *     to the cell size.
 */
TILEDB_EXPORT int32_t tiledb_attribute_set_fill_value(
    tiledb_ctx_t* ctx,
    tiledb_attribute_t* attr,
    const void* value,
    uint64_t size) TILEDB_NOEXCEPT;

/**
 * Gets the default fill value for the input attribute. This value will
 * be used for the input attribute whenever querying (1) an empty cell in
 * a dense array, or (2) a non-empty cell (in either dense or sparse array)
 * when values on the input attribute are missing (e.g., if the user writes
 * a subset of the attributes in a write operation).
 *
 * Applicable to both fixed-sized and var-sized attributes.
 *
 * **Example:**
 *
 * @code{.c}
 * // Assuming a int32 attribute
 * const int32_t* value;
 * uint64_t size;
 * tiledb_attribute_get_fill_value(ctx, attr, &value, &size);
 *
 * // Assuming a var char attribute
 * const char* value;
 * uint64_t size;
 * tiledb_attribute_get_fill_value(ctx, attr, &value, &size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The target attribute.
 * @param value A pointer to the fill value to get.
 * @param size The size of the fill value to get.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_get_fill_value(
    tiledb_ctx_t* ctx,
    tiledb_attribute_t* attr,
    const void** value,
    uint64_t* size) TILEDB_NOEXCEPT;

/**
 * Sets the default fill value for the input, nullable attribute. This value
 * will be used for the input attribute whenever querying (1) an empty cell in
 * a dense array, or (2) a non-empty cell (in either dense or sparse array)
 * when values on the input attribute are missing (e.g., if the user writes
 * a subset of the attributes in a write operation).
 *
 * Applicable to var-sized attributes.
 *
 * **Example:**
 *
 * @code{.c}
 * // Assumming a int32 attribute
 * int32_t value = 0;
 * uint64_t size = sizeof(value);
 * uint8_t valid = 0;
 * tiledb_attribute_set_fill_value_nullable(ctx, attr, &value, size, valid);
 *
 * // Assumming a var char attribute
 * const char* value = "foo";
 * uint64_t size = strlen(value);
 * uint8_t valid = 1;
 * tiledb_attribute_set_fill_value_nullable(ctx, attr, value, size, valid);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The target attribute.
 * @param value The fill value to set.
 * @param size The fill value size in bytes.
 * @param validity The validity fill value, zero for a null value and
 *     non-zero for a valid attribute.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note A call to `tiledb_attribute_cell_val_num` sets the fill value
 *     of the attribute to its default. Therefore, make sure you invoke
 *     `tiledb_attribute_set_fill_value_nullable` after deciding on the
 *     number of values this attribute will hold in each cell.
 *
 * @note For fixed-sized attributes, the input `size` should be equal
 *     to the cell size.
 */
TILEDB_EXPORT int32_t tiledb_attribute_set_fill_value_nullable(
    tiledb_ctx_t* ctx,
    tiledb_attribute_t* attr,
    const void* value,
    uint64_t size,
    uint8_t validity) TILEDB_NOEXCEPT;

/**
 * Gets the default fill value for the input, nullable attribute. This value
 * will be used for the input attribute whenever querying (1) an empty cell in
 * a dense array, or (2) a non-empty cell (in either dense or sparse array)
 * when values on the input attribute are missing (e.g., if the user writes
 * a subset of the attributes in a write operation).
 *
 * Applicable to both fixed-sized and var-sized attributes.
 *
 * **Example:**
 *
 * @code{.c}
 * // Assuming a int32 attribute
 * const int32_t* value;
 * uint64_t size;
 * uint8_t valid;
 * tiledb_attribute_get_fill_value_nullable(ctx, attr, &value, &size, &valid);
 *
 * // Assuming a var char attribute
 * const char* value;
 * uint64_t size;
 * uint8_t valid;
 * tiledb_attribute_get_fill_value_nullable(ctx, attr, &value, &size, &valid);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attr The target attribute.
 * @param value A pointer to the fill value to get.
 * @param size The size of the fill value to get.
 * @param valid The fill value validity to get.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_attribute_get_fill_value_nullable(
    tiledb_ctx_t* ctx,
    tiledb_attribute_t* attr,
    const void** value,
    uint64_t* size,
    uint8_t* valid) TILEDB_NOEXCEPT;

/* ********************************* */
/*               DOMAIN              */
/* ********************************* */

/**
 * Creates a TileDB domain.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_domain_t* domain;
 * tiledb_domain_alloc(ctx, &domain);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param domain The TileDB domain to be created.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_domain_alloc(
    tiledb_ctx_t* ctx, tiledb_domain_t** domain) TILEDB_NOEXCEPT;

/**
 * Destroys a TileDB domain, freeing associated memory.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_domain_t* domain;
 * tiledb_domain_alloc(ctx, &domain);
 * tiledb_domain_free(&domain);
 * @endcode
 *
 * @param domain The domain to be destroyed.
 */
TILEDB_EXPORT void tiledb_domain_free(tiledb_domain_t** domain) TILEDB_NOEXCEPT;

/**
 * Retrieves the domain's type.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_datatype_t type;
 * tiledb_domain_get_type(ctx, domain, &type);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param domain The domain.
 * @param type The type to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_domain_get_type(
    tiledb_ctx_t* ctx,
    const tiledb_domain_t* domain,
    tiledb_datatype_t* type) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of dimensions in a domain.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t dim_num;
 * tiledb_domain_get_ndim(ctx, domain, &dim_num);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param domain The domain
 * @param ndim The number of dimensions in a domain.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_domain_get_ndim(
    tiledb_ctx_t* ctx,
    const tiledb_domain_t* domain,
    uint32_t* ndim) TILEDB_NOEXCEPT;

/**
 * Adds a dimension to a TileDB domain.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_dimension_t* dim;
 * int64_t dim_domain[] = {1, 10};
 * int64_t tile_extent = 5;
 * tiledb_dimension_alloc(
 *     ctx, "dim_0", TILEDB_INT64, dim_domain, &tile_extent, &dim);
 * tiledb_domain_add_dimension(ctx, domain, dim);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param domain The domain to add the dimension to.
 * @param dim The dimension to be added.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_domain_add_dimension(
    tiledb_ctx_t* ctx,
    tiledb_domain_t* domain,
    tiledb_dimension_t* dim) TILEDB_NOEXCEPT;

/**
 * Retrieves a dimension object from a domain by index.
 *
 * **Example:**
 *
 * The following retrieves the first dimension from a domain.
 *
 * @code{.c}
 * tiledb_dimension_t* dim;
 * tiledb_domain_get_dimension_from_index(ctx, domain, 0, &dim);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param domain The domain to add the dimension to.
 * @param index The index of domain dimension
 * @param dim The retrieved dimension object.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_domain_get_dimension_from_index(
    tiledb_ctx_t* ctx,
    const tiledb_domain_t* domain,
    uint32_t index,
    tiledb_dimension_t** dim) TILEDB_NOEXCEPT;

/**
 * Retrieves a dimension object from a domain by name (key).
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_dimension_t* dim;
 * tiledb_domain_get_dimension_from_name(ctx, domain, "dim_0", &dim);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param domain The domain to add the dimension to.
 * @param name The name (key) of the requested dimension
 * @param dim The retrieved dimension object.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_domain_get_dimension_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_domain_t* domain,
    const char* name,
    tiledb_dimension_t** dim) TILEDB_NOEXCEPT;

/**
 * Checks whether the domain has a dimension of the given name.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t has_dim;
 * tiledb_domain_has_dimension(ctx, domain, "dim_0", &has_dim);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param domain The domain.
 * @param name The name of the dimension to check for.
 * @param has_dim Set to `1` if the domain has a dimension of the given name,
 *      else `0`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_domain_has_dimension(
    tiledb_ctx_t* ctx,
    const tiledb_domain_t* domain,
    const char* name,
    int32_t* has_dim) TILEDB_NOEXCEPT;

/**
 * Dumps the info of a domain in ASCII form to some output (e.g.,
 * file or `stdout`).
 *
 * **Example:**
 *
 * The following prints the domain dump to the standard output.
 *
 * @code{.c}
 * tiledb_domain_dump(ctx, domain, stdout);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param domain The domain.
 * @param out The output.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_domain_dump(
    tiledb_ctx_t* ctx,
    const tiledb_domain_t* domain,
    FILE* out) TILEDB_NOEXCEPT;

/* ********************************* */
/*             DIMENSION             */
/* ********************************* */

/**
 * Creates a dimension.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_dimension_t* dim;
 * int64_t dim_domain[] = {1, 10};
 * int64_t tile_extent = 5;
 * tiledb_dimension_alloc(
 *     ctx, "dim_0", TILEDB_INT64, dim_domain, &tile_extent, &dim);
 * @endcode
 *
 * Note: as laid out in the Storage Format,
 * the following Datatypes are not valid for Dimension:
 * TILEDB_CHAR, TILEDB_BLOB, TILEDB_BOOL, TILEDB_STRING_UTF8,
 * TILEDB_STRING_UTF16, TILEDB_STRING_UTF32, TILEDB_STRING_UCS2,
 * TILEDB_STRING_UCS4, TILEDB_ANY
 *
 * @param ctx The TileDB context.
 * @param name The dimension name.
 * @param type The dimension type.
 * @param dim_domain The dimension domain.
 * @param tile_extent The dimension tile extent.
 * @param dim The dimension to be created.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_alloc(
    tiledb_ctx_t* ctx,
    const char* name,
    tiledb_datatype_t type,
    const void* dim_domain,
    const void* tile_extent,
    tiledb_dimension_t** dim) TILEDB_NOEXCEPT;

/**
 * Destroys a TileDB dimension, freeing associated memory.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_dimension_free(&dim);
 * @endcode
 *
 * @param dim The dimension to be destroyed.
 */
TILEDB_EXPORT void tiledb_dimension_free(tiledb_dimension_t** dim)
    TILEDB_NOEXCEPT;

/**
 * Sets the filter list for a dimension.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_filter_list_alloc(ctx, &filter_list);
 * tiledb_filter_list_add_filter(ctx, filter_list, filter);
 * tiledb_dimension_set_filter_list(ctx, dim, filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param dim The target dimension.
 * @param filter_list The filter_list to be set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_set_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_dimension_t* dim,
    tiledb_filter_list_t* filter_list) TILEDB_NOEXCEPT;

/**
 * Sets the number of values per cell for a dimension. If this is not
 * used, the default is `1`.
 *
 * **Examples:**
 *
 * For a fixed-sized dimension:
 *
 * @code{.c}
 * tiledb_dimension_set_cell_val_num(ctx, dim, 3);
 * @endcode
 *
 * For a variable-sized dimension:
 *
 * @code{.c}
 * tiledb_dimension_set_cell_val_num(ctx, dim, TILEDB_VAR_NUM);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param dim The target dimension.
 * @param cell_val_num The number of values per cell.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_set_cell_val_num(
    tiledb_ctx_t* ctx,
    tiledb_dimension_t* dim,
    uint32_t cell_val_num) TILEDB_NOEXCEPT;

/**
 * Retrieves the filter list for a dimension.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_dimension_get_filter_list(ctx, dim, &filter_list);
 * tiledb_filter_list_free(&filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param dim The target dimension.
 * @param filter_list The filter list to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_get_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_dimension_t* dim,
    tiledb_filter_list_t** filter_list) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of values per cell for a dimension. For variable-sized
 * dimensions the result is TILEDB_VAR_NUM.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t num;
 * tiledb_dimension_get_cell_val_num(ctx, dim, &num);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param dim The dimension.
 * @param cell_val_num The number of values per cell to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_get_cell_val_num(
    tiledb_ctx_t* ctx,
    const tiledb_dimension_t* dim,
    uint32_t* cell_val_num) TILEDB_NOEXCEPT;

/**
 * Retrieves the dimension name.
 *
 * **Example:**
 *
 * @code{.c}
 * const char* dim_name;
 * tiledb_dimension_get_name(ctx, dim, &dim_name);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param dim The dimension.
 * @param name The name to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_get_name(
    tiledb_ctx_t* ctx,
    const tiledb_dimension_t* dim,
    const char** name) TILEDB_NOEXCEPT;

/**
 * Retrieves the dimension type.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_datatype_t dim_type;
 * tiledb_dimension_get_type(ctx, dim, &dim_type);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param dim The dimension.
 * @param type The type to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_get_type(
    tiledb_ctx_t* ctx,
    const tiledb_dimension_t* dim,
    tiledb_datatype_t* type) TILEDB_NOEXCEPT;

/**
 * Retrieves the domain of the dimension.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t* domain;
 * tiledb_dimension_get_domain(ctx, dim, &domain);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param dim The dimension.
 * @param domain The domain to be retrieved. Note that the defined type of
 *     input `domain` must be the same as the dimension type, otherwise the
 *     behavior is unpredictable (it will probably segfault).
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_get_domain(
    tiledb_ctx_t* ctx,
    const tiledb_dimension_t* dim,
    const void** domain) TILEDB_NOEXCEPT;

/**
 * Retrieves the tile extent of the dimension.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t* tile_extent;
 * tiledb_dimension_get_tile_extent(ctx, dim, &tile_extent);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param dim The dimension.
 * @param tile_extent The tile extent (pointer) to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_get_tile_extent(
    tiledb_ctx_t* ctx,
    const tiledb_dimension_t* dim,
    const void** tile_extent) TILEDB_NOEXCEPT;

/**
 * Dumps the contents of a dimension in ASCII form to some output (e.g.,
 * file or stdout).
 *
 * **Example:**
 *
 * The following prints the dimension dump to standard output.
 *
 * @code{.c}
 * tiledb_dimension_dump(ctx, dim, stdout);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param dim The dimension.
 * @param out The output.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_dimension_dump(
    tiledb_ctx_t* ctx,
    const tiledb_dimension_t* dim,
    FILE* out) TILEDB_NOEXCEPT;

/* ********************************* */
/*            ARRAY SCHEMA           */
/* ********************************* */

/**
 * Creates a TileDB array schema object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_t* array_schema;
 * tiledb_array_schema_alloc(ctx, TILEDB_DENSE, &array_schema);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_type The array type.
 * @param array_schema The TileDB array schema to be created.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_alloc(
    tiledb_ctx_t* ctx,
    tiledb_array_type_t array_type,
    tiledb_array_schema_t** array_schema) TILEDB_NOEXCEPT;

/**
 * Destroys an array schema, freeing associated memory.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_free(&array_schema);
 * @endcode
 *
 * @param array_schema The array schema to be destroyed.
 */
TILEDB_EXPORT void tiledb_array_schema_free(
    tiledb_array_schema_t** array_schema) TILEDB_NOEXCEPT;

/**
 * Adds an attribute to an array schema.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_attribute_t* attr;
 * tiledb_attribute_alloc(ctx, "my_attr", TILEDB_INT32, &attr);
 * tiledb_array_schema_add_attribute(ctx, array_schema, attr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param attr The attribute to be added.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_add_attribute(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_attribute_t* attr) TILEDB_NOEXCEPT;

/**
 * Sets whether the array can allow coordinate duplicates or not.
 * Applicable only to sparse arrays (it errors out if set to `1` for dense
 * arrays).
 *
 * **Example:**
 *
 * @code{.c}
 * int allows_dups = 1;
 * tiledb_array_schema_set_allows_dups(ctx, array_schema, allows_dups);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param allows_dups Whether or not the array allows coordinate duplicates.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_set_allows_dups(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    int allows_dups) TILEDB_NOEXCEPT;

/**
 * Gets whether the array can allow coordinate duplicates or not.
 * It should always be `0` for dense arrays.
 *
 * **Example:**
 *
 * @code{.c}
 * int allows_dups;
 * tiledb_array_schema_get_allows_dups(ctx, array_schema, &allows_dups);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param allows_dups Whether or not the array allows coordinate duplicates.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_allows_dups(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    int* allows_dups) TILEDB_NOEXCEPT;

/**
 * Returns the array schema version.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t version;
 * tiledb_array_schema_get_version(ctx, array_schema, &version);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param version The version.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_version(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    uint32_t* version) TILEDB_NOEXCEPT;

/**
 * Sets a domain for the array schema.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_domain_t* domain;
 * tiledb_domain_alloc(ctx, &domain);
 * // -- Add dimensions to the domain here -- //
 * tiledb_array_schema_set_domain(ctx, array_schema, domain);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param domain The domain to be set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_set_domain(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_domain_t* domain) TILEDB_NOEXCEPT;

/**
 * Sets the tile capacity.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_set_capacity(ctx, array_schema, 10000);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param capacity The capacity of a sparse data tile. Note that
 * sparse data tiles exist in sparse fragments, which can be created
 * in both sparse and dense arrays. For more details,
 * see [tutorials/tiling-sparse.html](tutorials/tiling-sparse.html).
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_set_capacity(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    uint64_t capacity) TILEDB_NOEXCEPT;

/**
 * Sets the cell order.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_set_cell_order(ctx, array_schema, TILEDB_ROW_MAJOR);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param cell_order The cell order to be set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_set_cell_order(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_layout_t cell_order) TILEDB_NOEXCEPT;

/**
 * Sets the tile order.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_set_cell_order(ctx, array_schema, TILEDB_COL_MAJOR);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param tile_order The tile order to be set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_set_tile_order(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_layout_t tile_order) TILEDB_NOEXCEPT;

/**
 * Sets the filter list to use for the coordinates.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_filter_list_alloc(ctx, &filter_list);
 * tiledb_filter_list_add_filter(ctx, filter_list, filter);
 * tiledb_array_schema_set_coords_filter_list(ctx, array_schema, filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param filter_list The filter list to be set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_set_coords_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_filter_list_t* filter_list) TILEDB_NOEXCEPT;

/**
 * Sets the filter list to use for the offsets of variable-sized attribute
 * values.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_filter_list_alloc(ctx, &filter_list);
 * tiledb_filter_list_add_filter(ctx, filter_list, filter);
 * tiledb_array_schema_set_offsets_filter_list(ctx, array_schema, filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param filter_list The filter list to be set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_set_offsets_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_filter_list_t* filter_list) TILEDB_NOEXCEPT;

/**
 * Sets the filter list to use for the validity array of nullable attribute
 * values.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_filter_list_alloc(ctx, &filter_list);
 * tiledb_filter_list_add_filter(ctx, filter_list, filter);
 * tiledb_array_schema_set_validity_filter_list(ctx, array_schema, filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param filter_list The filter list to be set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_set_validity_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_filter_list_t* filter_list) TILEDB_NOEXCEPT;

/**
 * Checks the correctness of the array schema.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_check(ctx, array_schema);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @return `TILEDB_OK` if the array schema is correct and `TILEDB_ERR` upon any
 *     error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_check(
    tiledb_ctx_t* ctx, tiledb_array_schema_t* array_schema) TILEDB_NOEXCEPT;

/**
 * Retrieves the schema of an array from the disk, creating an array schema
 * struct.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_t* array_schema;
 * tiledb_array_schema_load(ctx, "s3://tiledb_bucket/my_array", &array_schema);
 * // Make sure to free the array schema in the end
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The array whose schema will be retrieved.
 * @param array_schema The array schema to be retrieved, or `NULL` upon error.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_load(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_array_schema_t** array_schema) TILEDB_NOEXCEPT;

/**
 * Retrieves the schema of an encrypted array from the disk, creating an array
 * schema struct.
 *
 * **Example:**
 *
 * @code{.c}
 * // Load AES-256 key from disk, environment variable, etc.
 * uint8_t key[32] = ...;
 * tiledb_array_schema_t* array_schema;
 * tiledb_array_schema_load_with_key(
 *     ctx, "s3://tiledb_bucket/my_array", TILEDB_AES_256_GCM,
 *     key, sizeof(key), &array_schema);
 * // Make sure to free the array schema in the end
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The array whose schema will be retrieved.
 * @param encryption_type The encryption type to use.
 * @param encryption_key The encryption key to use.
 * @param key_length Length in bytes of the encryption key.
 * @param array_schema The array schema to be retrieved, or `NULL` upon error.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_schema_load_with_key(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_encryption_type_t encryption_type,
    const void* encryption_key,
    uint32_t key_length,
    tiledb_array_schema_t** array_schema) TILEDB_NOEXCEPT;

/**
 * Retrieves the array type.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_t* array_schema;
 * tiledb_array_schema_load(ctx, "s3://tiledb_bucket/my_array", array_schema);
 * tiledb_array_type_t* array_type;
 * tiledb_array_schema_get_array_type(ctx, array_schema, &array_type);
 * // Make sure to free the array schema in the end
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param array_type The array type to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_array_type(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    tiledb_array_type_t* array_type) TILEDB_NOEXCEPT;

/**
 * Retrieves the capacity.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t capacity;
 * tiledb_array_schema_get_capacity(ctx, array_schema, &capacity);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param capacity The capacity to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_capacity(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    uint64_t* capacity) TILEDB_NOEXCEPT;

/**
 * Retrieves the cell order.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_layout_t cell_order;
 * tiledb_array_schema_get_cell_order(ctx, array_schema, &cell_order);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param cell_order The cell order to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_cell_order(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    tiledb_layout_t* cell_order) TILEDB_NOEXCEPT;

/**
 * Retrieves the filter list used for the coordinates.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_array_schema_get_coords_filter_list(ctx, array_schema, &filter_list);
 * tiledb_filter_list_free(ctx, &filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param filter_list The filter list to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_coords_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_filter_list_t** filter_list) TILEDB_NOEXCEPT;

/**
 * Retrieves the filter list used for the offsets.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_array_schema_get_offsets_filter_list(ctx, array_schema, &filter_list);
 * tiledb_filter_list_free(ctx, &filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param filter_list The filter list to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_offsets_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_filter_list_t** filter_list) TILEDB_NOEXCEPT;

/**
 * Retrieves the filter list used for validity maps.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_filter_list_t* filter_list;
 * tiledb_array_schema_get_validity_filter_list(ctx, array_schema,
 * &filter_list); tiledb_filter_list_free(ctx, &filter_list);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param filter_list The filter list to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_validity_filter_list(
    tiledb_ctx_t* ctx,
    tiledb_array_schema_t* array_schema,
    tiledb_filter_list_t** filter_list) TILEDB_NOEXCEPT;

/**
 * Retrieves the array domain.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_domain_t* domain;
 * tiledb_array_schema_get_domain(ctx, array_schema, &domain);
 * // Make sure to delete domain in the end
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param domain The array domain to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_domain(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    tiledb_domain_t** domain) TILEDB_NOEXCEPT;

/**
 * Retrieves the tile order.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_layout_t tile_order;
 * tiledb_array_schema_get_tile_order(ctx, array_schema, &tile_order);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param tile_order The tile order to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_tile_order(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    tiledb_layout_t* tile_order) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of array attributes.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t attr_num;
 * tiledb_array_schema_get_attribute_num(ctx, array_schema, &attr_num);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param attribute_num The number of attributes to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_attribute_num(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    uint32_t* attribute_num) TILEDB_NOEXCEPT;

/**
 * Retrieves an attribute given its index.
 *
 * Attributes are ordered the same way they were defined
 * when constructing the array schema.
 *
 * **Example:**
 *
 * The following retrieves the first attribute in the schema.
 *
 * @code{.c}
 * tiledb_attribute_t* attr;
 * tiledb_array_schema_get_attribute_from_index(ctx, array_schema, 0, &attr);
 * // Make sure to delete the retrieved attribute in the end.
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param index The index of the attribute to retrieve.
 * @param attr The attribute object to retrieve.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_attribute_from_index(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    uint32_t index,
    tiledb_attribute_t** attr) TILEDB_NOEXCEPT;

/**
 * Retrieves an attribute given its name (key).
 *
 * **Example:**
 *
 * The following retrieves the first attribute in the schema.
 *
 * @code{.c}
 * tiledb_attribute_t* attr;
 * tiledb_array_schema_get_attribute_from_name(
 *     ctx, array_schema, "attr_0", &attr);
 * // Make sure to delete the retrieved attribute in the end.
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param name The name (key) of the attribute to retrieve.
 * @param attr THe attribute object to retrieve.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_get_attribute_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    const char* name,
    tiledb_attribute_t** attr) TILEDB_NOEXCEPT;

/**
 * Checks whether the array schema has an attribute of the given name.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t has_attr;
 * tiledb_array_schema_has_attribute(ctx, array_schema, "attr_0", &has_attr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param name The name of the attribute to check for.
 * @param has_attr Set to `1` if the array schema has an attribute of the
 *      given name, else `0`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_has_attribute(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    const char* name,
    int32_t* has_attr) TILEDB_NOEXCEPT;

/**
 * Dumps the array schema in ASCII format in the selected output.
 *
 * **Example:**
 *
 * The following prints the array schema dump in standard output.
 *
 * @code{.c}
 * tiledb_array_schema_dump(ctx, array_schema, stdout);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_schema The array schema.
 * @param out The output.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_schema_dump(
    tiledb_ctx_t* ctx,
    const tiledb_array_schema_t* array_schema,
    FILE* out) TILEDB_NOEXCEPT;

/* ********************************* */
/*               QUERY               */
/* ********************************* */

/**
 * Creates a TileDB query object. Note that the query object is associated
 * with a specific array object. The query type (read or write) is inferred
 * from the array object, which was opened with a specific query type.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "file:///my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_WRITE);
 * tiledb_query_t* query;
 * tiledb_query_alloc(ctx, array, TILEDB_WRITE, &query);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query object to be created.
 * @param array An opened array object.
 * @param query_type The query type. This must comply with the query type
 *     `array` was opened.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_alloc(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    tiledb_query_type_t query_type,
    tiledb_query_t** query) TILEDB_NOEXCEPT;

/**
 * Retrieves the stats from a Query.
 *
 * **Example:**
 *
 * @code{.c}
 * char* stats_json;
 * tiledb_query_get_stats(ctx, query, &stats_json);
 * // Make sure to free the retrieved `stats_json`
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query object.
 * @param stats_json The output json. The caller takes ownership
 *   of the c-string.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_stats(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    char** stats_json) TILEDB_NOEXCEPT;

/**
 * Set the query config
 *
 * Setting the query config will also set the subarray configuration in order to
 * maintain existing behavior. If you wish the subarray to have a different
 * configuration than the query, set it after calling tiledb_query_set_config.
 *
 * Setting the configuration with this function overrides the following
 * Query-level parameters only:
 *
 * - `sm.memory_budget`
 * - `sm.memory_budget_var`
 * - `sm.var_offsets.mode`
 * - `sm.var_offsets.extra_element`
 * - `sm.var_offsets.bitsize`
 * - `sm.check_coord_dups`
 * - `sm.check_coord_oob`
 * - `sm.check_global_order`
 * - `sm.dedup_coords`
 */
TILEDB_EXPORT int32_t tiledb_query_set_config(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    tiledb_config_t* config) TILEDB_NOEXCEPT;

/**
 * Retrieves the config from a Query.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_config_t* config;
 * tiledb_query_get_config(ctx, vfs, &config);
 * // Make sure to free the retrieved config
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query object.
 * @param config The config to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_config(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    tiledb_config_t** config) TILEDB_NOEXCEPT;
/**
 * Indicates that the query will write or read a subarray, and provides
 * the appropriate information.
 *
 * **Example:**
 *
 * The following sets a 2D subarray [0,10], [20, 30] to the query.
 *
 * @code{.c}
 * uint64_t subarray[] = { 0, 10, 20, 30};
 * tiledb_query_set_subarray(ctx, query, subarray);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param subarray The subarray in which the array read/write will be
 *     constrained on. It should be a sequence of [low, high] pairs (one
 *     pair per dimension). For the case of writes, this is meaningful only
 *     for dense arrays, and specifically dense writes. Note that `subarray`
 *     must have the same type as the domain.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 *
 * @note If you set the subarray of a completed, incomplete or in-progress
 *     query, this function will clear the internal state and render it
 *     as uninitialized. However, the potentially set layout and attribute
 *     buffers will be retained. This is useful when the user wishes to
 *     fix the attributes and layout, but explore different subarrays with
 *     the same `tiledb_query_t` object (i.e., without having to create
 *     a new object).
 *
 * @note This function will error in the following case, provided that
 *     this is a write query:
 *     (i) the array is dense and the
 *     layout has been set to `TILEDB_UNORDERED`. In this case,
 *     if the user sets the layout to `TILEDB_UNORDERED` **after**
 *     the subarray has been set, the subarray will simply be ignored.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_set_subarray(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const void* subarray) TILEDB_NOEXCEPT;

/**
 * Indicates that the query will write or read a subarray, and provides
 * the appropriate information.
 *
 * **Example:**
 *
 * The following sets a 2D subarray [0,10], [20, 30] to the query.
 *
 * @code{.c}
 * tiledb_subarray_t *subarray;
 * tiledb_subarray_alloc(ctx, array, &subarray);
 * uint64_t subarray_v[] = { 0, 10, 20, 30};
 * tiledb_subarray_set_subarray(ctx, subarray, subarray_v);
 * tiledb_query_set_subarray_t(ctx, query, subarray);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param subarray The subarray by which the array read/write will be
 *     constrained. It should be a sequence of [low, high] pairs (one
 *     pair per dimension). For the case of writes, this is meaningful only
 *     for dense arrays, and specifically dense writes. Note that `subarray`
 *     must have the same type as the domain.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 *
 * @note If you set the subarray of a completed, or uninitialized
 *     query, this function will clear the internal state and render it
 *     as uninitialized. However, the potentially set layout and attribute
 *     buffers will be retained. This is useful when the user wishes to
 *     fix the attributes and layout, but explore different subarrays with
 *     the same `tiledb_query_t` object (i.e., without having to create
 *     a new object).
 *
 * @note Setting the subarray in sparse writes is meaningless and, thus,
 *     this function will error in the following two cases, provided that
 *     this is a write query:
 *     (i) the array is sparse, and (ii) the array is dense and the
 *     layout has been set to `TILEDB_UNORDERED`. In the second case,
 *     if the user sets the layout to `TILEDB_UNORDERED` **after**
 *     the subarray has been set, the subarray will simply be ignored.
 */
TILEDB_EXPORT int32_t tiledb_query_set_subarray_t(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const tiledb_subarray_t* subarray) TILEDB_NOEXCEPT;

/**
 * Sets the buffer for a fixed-sized attribute/dimension to a query, which will
 * either hold the values to be written (if it is a write query), or will hold
 * the results from a read query.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t a1[100];
 * uint64_t a1_size = sizeof(a1);
 * tiledb_query_set_buffer(ctx, query, "a1", a1, &a1_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to set the buffer for. Note that
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer The buffer that either have the input data to be written,
 *     or will hold the data to be read.
 * @param buffer_size In the case of writes, this is the size of `buffer`
 *     in bytes. In the case of reads, this initially contains the allocated
 *     size of `buffer`, but after the termination of the query
 *     it will contain the size of the useful (read) data in `buffer`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_set_buffer(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    void* buffer,
    uint64_t* buffer_size) TILEDB_NOEXCEPT;

/**
 * Sets the buffer for a var-sized attribute/dimension to a query, which will
 * either hold the values to be written (if it is a write query), or will hold
 * the results from a read query.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t a2_off[10];
 * uint64_t a2_off_size = sizeof(a2_off);
 * char a2_val[100];
 * uint64_t a2_val_size = sizeof(a2_val);
 * tiledb_query_set_buffer_var(
 *     ctx, query, "a2", a2_off, &a2_off_size, a2_val, &a2_val_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to set the buffer for.
 * @param buffer_off The buffer that either have the input data to be written,
 *     or will hold the data to be read. This buffer holds the starting offsets
 *     of each cell value in `buffer_val`.
 * @param buffer_off_size In the case of writes, it is the size of `buffer_off`
 *     in bytes. In the case of reads, this initially contains the allocated
 *     size of `buffer_off`, but after the *end of the query*
 *     (`tiledb_query_submit`) it will contain the size of the useful (read)
 *     data in `buffer_off`.
 * @param buffer_val The buffer that either have the input data to be written,
 *     or will hold the data to be read. This buffer holds the actual var-sized
 *     cell values.
 * @param buffer_val_size In the case of writes, it is the size of `buffer_val`
 *     in bytes. In the case of reads, this initially contains the allocated
 *     size of `buffer_val`, but after the termination of the function
 *     it will contain the size of the useful (read) data in `buffer_val`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_set_buffer_var(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    uint64_t* buffer_off,
    uint64_t* buffer_off_size,
    void* buffer_val,
    uint64_t* buffer_val_size) TILEDB_NOEXCEPT;

/**
 * Sets the buffer for a fixed-sized, nullable attribute to a query, which will
 * either hold the values to be written (if it is a write query), or will hold
 * the results from a read query. The validity buffer is a byte map, where each
 * non-zero byte represents a valid (i.e. "non-null") attribute value.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t a1[100];
 * uint64_t a1_size = sizeof(a1);
 * uint8_t a1_validity[100];
 * uint64_t a1_validity_size = sizeof(a1_validity);
 * tiledb_query_set_buffer_nullable(
 *   ctx, query, "a1", a1, &a1_size, a1_validity, &a1_validity_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to set the buffer for. Note that
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer The buffer that either have the input data to be written,
 *     or will hold the data to be read.
 * @param buffer_size In the case of writes, this is the size of `buffer`
 *     in bytes. In the case of reads, this initially contains the allocated
 *     size of `buffer`, but after the termination of the query
 *     it will contain the size of the useful (read) data in `buffer`.
 * @param buffer_validity_bytemap The validity byte map that has exactly
 *     one value for each value in `buffer`.
 * @param buffer_validity_bytemap_size In the case of writes, this is the
 *     size of `buffer_validity_bytemap` in bytes. In the case of reads,
 *     this initially contains the allocated size of `buffer_validity_bytemap`,
 *     but after the termination of the query it will contain the size of the
 *     useful (read) data in `buffer_validity_bytemap`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_set_buffer_nullable(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    void* buffer,
    uint64_t* buffer_size,
    uint8_t* buffer_validity_bytemap,
    uint64_t* buffer_validity_bytemap_size) TILEDB_NOEXCEPT;

/**
 * Sets the buffer for a var-sized, nullable attribute to a query, which will
 * either hold the values to be written (if it is a write query), or will hold
 * the results from a read query.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t a2_off[10];
 * uint64_t a2_off_size = sizeof(a2_off);
 * char a2_val[100];
 * uint64_t a2_val_size = sizeof(a2_val);
 * uint8_t a2_validity[100];
 * uint64_t a2_validity_size = sizeof(a2_validity);
 * tiledb_query_set_buffer_var(
 *     ctx, query, "a2", a2_off, &a2_off_size, a2_val, &a2_val_size,
 *     a2_validity, &a2_validity_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to set the buffer for.
 * @param buffer_off The buffer that either have the input data to be written,
 *     or will hold the data to be read. This buffer holds the starting offsets
 *     of each cell value in `buffer_val`.
 * @param buffer_off_size In the case of writes, it is the size of `buffer_off`
 *     in bytes. In the case of reads, this initially contains the allocated
 *     size of `buffer_off`, but after the *end of the query*
 *     (`tiledb_query_submit`) it will contain the size of the useful (read)
 *     data in `buffer_off`.
 * @param buffer_val The buffer that either have the input data to be written,
 *     or will hold the data to be read. This buffer holds the actual var-sized
 *     cell values.
 * @param buffer_val_size In the case of writes, it is the size of `buffer_val`
 *     in bytes. In the case of reads, this initially contains the allocated
 *     size of `buffer_val`, but after the termination of the function
 *     it will contain the size of the useful (read) data in `buffer_val`.
 * @param buffer_validity_bytemap The validity byte map that has exactly
 *     one value for each value in `buffer`.
 * @param buffer_validity_bytemap_size In the case of writes, this is the
 *     size of `buffer_validity_bytemap` in bytes. In the case of reads,
 *     this initially contains the allocated size of `buffer_validity_bytemap`,
 *     but after the termination of the query it will contain the size of the
 *     useful (read) data in `buffer_validity_bytemap`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_set_buffer_var_nullable(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    uint64_t* buffer_off,
    uint64_t* buffer_off_size,
    void* buffer_val,
    uint64_t* buffer_val_size,
    uint8_t* buffer_validity_bytemap,
    uint64_t* buffer_validity_bytemap_size) TILEDB_NOEXCEPT;

/**
 * Sets the buffer for an attribute/dimension to a query, which will
 * either hold the values to be written (if it is a write query), or will hold
 * the results from a read query.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t a1[100];
 * uint64_t a1_size = sizeof(a1);
 * tiledb_query_set_data_buffer(ctx, query, "a1", a1, &a1_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to set the buffer for. Note that
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer The buffer that either have the input data to be written,
 *     or will hold the data to be read.
 * @param buffer_size In the case of writes, this is the size of `buffer`
 *     in bytes. In the case of reads, this initially contains the allocated
 *     size of `buffer`, but after the termination of the query
 *     it will contain the size of the useful (read) data in `buffer`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_set_data_buffer(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    void* buffer,
    uint64_t* buffer_size) TILEDB_NOEXCEPT;

/**
 * Sets the starting offsets of each cell value in the data buffer.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t a1[100];
 * uint64_t a1_size = sizeof(a1);
 * tiledb_query_set_offsets_buffer(ctx, query, "a1", a1, &a1_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to set the buffer for. Note that
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer This buffer holds the starting offsets
 *     of each cell value in `buffer_val`.
 * @param buffer_size In the case of writes, it is the size of `buffer_off`
 *     in bytes. In the case of reads, this initially contains the allocated
 *     size of `buffer_off`, but after the *end of the query*
 *     (`tiledb_query_submit`) it will contain the size of the useful (read)
 *     data in `buffer_off`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_set_offsets_buffer(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    uint64_t* buffer,
    uint64_t* buffer_size) TILEDB_NOEXCEPT;

/**
 * Sets the validity byte map that has exactly one value for each value in the
 * data buffer.
 *
 * **Example:**
 *
 * @code{.c}
 * uint8_t a1[100];
 * uint64_t a1_size = sizeof(a1);
 * tiledb_query_set_validity_buffer(ctx, query, "a1", a1, &a1_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to set the buffer for. Note that
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer The validity byte map that has exactly
 *     one value for each value in `buffer`.
 * @param buffer_size In the case of writes, this is the
 *     size of `buffer_validity_bytemap` in bytes. In the case of reads,
 *     this initially contains the allocated size of `buffer_validity_bytemap`,
 *     but after the termination of the query it will contain the size of the
 *     useful (read) data in `buffer_validity_bytemap`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_set_validity_buffer(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    uint8_t* buffer,
    uint64_t* buffer_size) TILEDB_NOEXCEPT;

/**
 * Gets the buffer of a fixed-sized attribute/dimension from a query. If the
 * buffer has not been set, then `buffer` is set to `nullptr`.
 *
 * **Example:**
 *
 * @code{.c}
 * int* a1;
 * uint64_t* a1_size;
 * tiledb_query_get_buffer(ctx, query, "a1", &a1, &a1_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to get the buffer for. Note that the
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer The buffer to retrieve.
 * @param buffer_size A pointer to the size of the buffer. Note that this is
 *     a double pointer and returns the original variable address from
 *     `set_buffer`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_buffer(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    void** buffer,
    uint64_t** buffer_size) TILEDB_NOEXCEPT;

/**
 * Gets the values and offsets buffers for a var-sized attribute/dimension
 * to a query. If the buffers have not been set, then `buffer_off` and
 * `buffer_val` are set to `nullptr`.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t* a2_off;
 * uint64_t* a2_off_size;
 * char* a2_val;
 * uint64_t* a2_val_size;
 * tiledb_query_get_buffer_var(
 *     ctx, query, "a2", &a2_off, &a2_off_size, &a2_val, &a2_val_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to set the buffer for.
 * @param buffer_off The offsets buffer to be retrieved.
 * @param buffer_off_size A pointer to the size of the offsets buffer. Note that
 *     this is a `uint_64**` pointer and returns the original variable address
 * from `set_buffer`.
 * @param buffer_val The values buffer to be retrieved.
 * @param buffer_val_size A pointer to the size of the values buffer. Note that
 *     this is a `uint_64**` pointer and returns the original variable address
 * from `set_buffer`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_buffer_var(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    uint64_t** buffer_off,
    uint64_t** buffer_off_size,
    void** buffer_val,
    uint64_t** buffer_val_size) TILEDB_NOEXCEPT;

/**
 * Gets the buffer of a fixed-sized, nullable attribute from a query. If the
 * buffer has not been set, then `buffer` and `buffer_validity_bytemap` are
 * set to `nullptr`.
 *
 * **Example:**
 *
 * @code{.c}
 * int* a1;
 * uint64_t* a1_size;
 * uint8_t* a1_validity;
 * uint64_t* a1_validity_size;
 * tiledb_query_get_buffer_nullable(
 *   ctx, query, "a1", &a1, &a1_size, &a1_validity, &a1_validity_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to get the buffer for. Note that the
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer The buffer to retrieve.
 * @param buffer_size A pointer to the size of the buffer. Note that this is
 *     a double pointer and returns the original variable address from
 *     `set_buffer`.
 * @param buffer_validity_bytemap The validity bytemap buffer to retrieve.
 * @param buffer_validity_bytemap_size A pointer to the size of the validity
 *     bytemap buffer. Note that this is a double pointer and returns the
 * origina variable address from `set_buffer_nullable`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_buffer_nullable(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    void** buffer,
    uint64_t** buffer_size,
    uint8_t** buffer_validity_bytemap,
    uint64_t** buffer_validity_bytemap_size) TILEDB_NOEXCEPT;

/**
 * Gets the values and offsets buffers for a var-sized, nullable attribute
 * to a query. If the buffers have not been set, then `buffer_off`,
 * `buffer_val`, and `buffer_validity_bytemap` are set to `nullptr`.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t* a2_off;
 * uint64_t* a2_off_size;
 * char* a2_val;
 * uint64_t* a2_val_size;
 * uint8_t* a2_validity;
 * uint64_t* a2_validity_size;
 * tiledb_query_get_buffer_var_nullable(
 *     ctx, query, "a2", &a2_off, &a2_off_size, &a2_val, &a2_val_size,
 *     &a2_validity, &a2_validity_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to set the buffer for.
 * @param buffer_off The offsets buffer to be retrieved.
 * @param buffer_off_size A pointer to the size of the offsets buffer. Note that
 *     this is a `uint_64**` pointer and returns the original variable address
 * from `set_buffer`.
 * @param buffer_val The values buffer to be retrieved.
 * @param buffer_val_size A pointer to the size of the values buffer. Note that
 *     this is a `uint_64**` pointer and returns the original variable address
 * from `set_buffer`.
 * @param buffer_validity_bytemap The validity bytemap buffer to retrieve.
 * @param buffer_validity_bytemap_size A pointer to the size of the validity
 *     bytemap buffer. Note that this is a double pointer and returns the
 * origina variable address from `set_buffer_var_nullable`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_buffer_var_nullable(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    uint64_t** buffer_off,
    uint64_t** buffer_off_size,
    void** buffer_val,
    uint64_t** buffer_val_size,
    uint8_t** buffer_validity_bytemap,
    uint64_t** buffer_validity_bytemap_size) TILEDB_NOEXCEPT;

/**
 * Gets the buffer of a fixed-sized attribute/dimension from a query. If the
 * buffer has not been set, then `buffer` is set to `nullptr`.
 *
 * **Example:**
 *
 * @code{.c}
 * int* a1;
 * uint64_t* a1_size;
 * tiledb_query_get_data_buffer(ctx, query, "a1", &a1, &a1_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to get the buffer for. Note that the
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer The buffer to retrieve.
 * @param buffer_size A pointer to the size of the buffer. Note that this is
 *     a double pointer and returns the original variable address from
 *     `set_buffer`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_data_buffer(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    void** buffer,
    uint64_t** buffer_size) TILEDB_NOEXCEPT;

/**
 * Gets the starting offsets of each cell value in the data buffer.
 *
 * **Example:**
 *
 * @code{.c}
 * int* a1;
 * uint64_t* a1_size;
 * tiledb_query_get_offsets_buffer(ctx, query, "a1", &a1, &a1_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to get the buffer for. Note that the
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer The buffer to retrieve.
 * @param buffer_size A pointer to the size of the buffer. Note that this is
 *     a double pointer and returns the original variable address from
 *     `set_buffer`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_offsets_buffer(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    uint64_t** buffer,
    uint64_t** buffer_size) TILEDB_NOEXCEPT;

/**
 * Gets the validity byte map that has exactly one value for each value in the
 * data buffer.
 *
 * **Example:**
 *
 * @code{.c}
 * int* a1;
 * uint64_t* a1_size;
 * tiledb_query_get_validity_buffer(ctx, query, "a1", &a1, &a1_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param name The attribute/dimension to get the buffer for. Note that the
 *     zipped coordinates have special name `TILEDB_COORDS`.
 * @param buffer The buffer to retrieve.
 * @param buffer_size A pointer to the size of the buffer. Note that this is
 *     a double pointer and returns the original variable address from
 *     `set_buffer`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_validity_buffer(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* name,
    uint8_t** buffer,
    uint64_t** buffer_size) TILEDB_NOEXCEPT;

/**
 * Sets the layout of the cells to be written or read.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_set_layout(ctx, query, TILEDB_ROW_MAJOR);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param layout For a write query, this specifies the order of the cells
 *     provided by the user in the buffers. For a read query, this specifies
 *     the order of the cells that will be retrieved as results and stored
 *     in the user buffers. The layout can be one of the following:
 *    - `TILEDB_COL_MAJOR`:
 *      This means column-major order with respect to the subarray.
 *    - `TILEDB_ROW_MAJOR`:
 *      This means row-major order with respect to the subarray.
 *    - `TILEDB_GLOBAL_ORDER`:
 *      This means that cells are stored or retrieved in the array global
 *      cell order.
 *    - `TILEDB_UNORDERED`:
 *      This is applicable only to reads and writes for sparse arrays, or for
 *      sparse writes to dense arrays. For writes, it specifies that the cells
 *      are unordered and, hence, TileDB must sort the cells in the global cell
 *      order prior to writing. For reads, TileDB will return the cells without
 *      any particular order, which will often lead to better performance.
 * * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_set_layout(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    tiledb_layout_t layout) TILEDB_NOEXCEPT;

/**
 * Sets the query condition to be applied on a read.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_condition_t* query_condition;
 * tiledb_query_condition_alloc(ctx, &query_condition);
 * uint32_t value = 5;
 * tiledb_query_condition_init(
 *   ctx, query_condition, "longitude", &value, sizeof(value), TILEDB_LT);
 * tiledb_query_set_condition(ctx, query, query_condition);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The TileDB query.
 * @param cond The TileDB query condition.
 */
TILEDB_EXPORT int32_t tiledb_query_set_condition(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const tiledb_query_condition_t* cond) TILEDB_NOEXCEPT;

/**
 * Flushes all internal state of a query object and finalizes the query.
 * This is applicable only to global layout writes. It has no effect for
 * any other query type.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_t* query;
 * // ... Your code here ... //
 * tiledb_query_finalize(ctx, query);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query object to be flushed.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t
tiledb_query_finalize(tiledb_ctx_t* ctx, tiledb_query_t* query) TILEDB_NOEXCEPT;

/**
 * Submits and finalizes the query.
 * This is applicable only to global layout writes. The function will
 * error out if called on a query with non global layout.
 * Its purpose is to submit the final chunk (partial or full tile) in
 * a global order write query.
 * `tiledb_query_submit_and_finalize` drops the tile alignment restriction
 * of the buffers (i.e. compared to the regular global layout submit call)
 * given the last chunk of a global order write is most frequently smaller
 * in size than a tile.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_t* query;
 * while (stop_condition) {
 *   tiledb_query_set_buffer(ctx, query, attr, tile_aligned_buffer, &size);
 *   tiledb_query_submit(ctx, query);
 * }
 * tiledb_query_set_buffer(ctx, query, attr, final_chunk, &size);
 * tiledb_query_submit_and_finalize(ctx, query);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query object to be flushed.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_submit_and_finalize(
    tiledb_ctx_t* ctx, tiledb_query_t* query) TILEDB_NOEXCEPT;

/**
 * Frees a TileDB query object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_free(&query);
 * @endcode
 *
 * @param query The query object to be deleted.
 */
TILEDB_EXPORT void tiledb_query_free(tiledb_query_t** query) TILEDB_NOEXCEPT;

/**
 * Submits a TileDB query.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_submit(ctx, query);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query to be submitted.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note `tiledb_query_finalize` must be invoked after finish writing in
 *     global layout (via repeated invocations of `tiledb_query_submit`),
 *     in order to flush any internal state.
 *
 * @note For the case of reads, if the returned status is `TILEDB_INCOMPLETE`,
 *    TileDB could not fit the entire result in the user's buffers. In this
 *    case, the user should consume the read results (if any), optionally
 *    reset the buffers with `tiledb_query_set_buffer`, and then resubmit the
 *    query until the status becomes `TILEDB_COMPLETED`. If all buffer sizes
 *    after the termination of this function become 0, then this means that
 *    **no** useful data was read into the buffers, implying that larger
 *    buffers are needed for the query to proceed. In this case, the users
 *    must reallocate their buffers (increasing their size), reset the buffers
 *    with `tiledb_query_set_buffer`, and resubmit the query.
 */
TILEDB_EXPORT int32_t
tiledb_query_submit(tiledb_ctx_t* ctx, tiledb_query_t* query) TILEDB_NOEXCEPT;

/**
 * Submits a TileDB query in asynchronous mode.
 *
 * **Examples:**
 *
 * Submit without a callback.
 *
 * @code{.c}
 * tiledb_query_submit_async(ctx, query, NULL, NULL);
 * @endcode
 *
 * Submit with a callback function `print` that takes as input message
 * `msg` and prints it upon completion of the query.
 *
 * @code{.c}
 * const char* msg = "Query completed";
 * tiledb_query_submit_async(ctx, &query, foo, msg);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query to be submitted.
 * @param callback The function to be called when the query completes.
 * @param callback_data The data to be passed to the \p callback function.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 *
 * @note `tiledb_query_finalize` must be invoked after finish writing in
 *     global layout (via repeated invocations of `tiledb_query_submit`),
 *     in order to flush any internal state.
 *
 * @note For the case of reads, if the returned status is `TILEDB_INCOMPLETE`,
 *    TileDB could not fit the entire result in the user's buffers. In this
 *    case, the user should consume the read results (if any), optionally
 *    reset the buffers with `tiledb_query_set_buffer`, and then resubmit the
 *    query until the status becomes `TILEDB_COMPLETED`. If all buffer sizes
 *    after the termination of this function become 0, then this means that
 *    **no** useful data was read into the buffers, implying that larger
 *    buffers are needed for the query to proceed. In this case, the users
 *    must reallocate their buffers (increasing their size), reset the buffers
 *    with `tiledb_query_set_buffer`, and resubmit the query.
 *
 * @note \p callback will be executed in a thread managed by TileDB's internal
 *    thread pool. To allow TileDB to reuse the thread and avoid starving the
 *    thread pool, long-running callbacks should be dispatched to another
 *    thread.
 */
TILEDB_EXPORT int32_t tiledb_query_submit_async(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    void (*callback)(void*),
    void* callback_data) TILEDB_NOEXCEPT;

/**
 * Checks if the query has returned any results. Applicable only to
 * read queries; it sets `has_results` to `0 in the case of writes.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t has_results;
 * tiledb_query_has_results(ctx, query, &has_results);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query.
 * @param has_results Set to `1` if the query returned results and `0`
 *     otherwise.
 * @return `TILEDB_OK` upon success, and `TILEDB_ERR` upon error.
 */
TILEDB_EXPORT int32_t tiledb_query_has_results(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    int32_t* has_results) TILEDB_NOEXCEPT;

/**
 * Retrieves the status of a query.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_status_t status;
 * tiledb_query_get_status(ctx, query, &status);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query.
 * @param status The query status to be retrieved.
 * @return `TILEDB_OK` upon success, and `TILEDB_ERR` upon error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_status(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    tiledb_query_status_t* status) TILEDB_NOEXCEPT;

/**
 * Retrieves the query type.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_type_t query_type;
 * tiledb_query_get_status(ctx, query, &query_type);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query.
 * @param query_type The query type to be retrieved.
 * @return `TILEDB_OK` upon success, and `TILEDB_ERR` upon error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_type(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    tiledb_query_type_t* query_type) TILEDB_NOEXCEPT;

/**
 * Retrieves the query layout.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_layout_t query_layout;
 * tiledb_query_get_layout(ctx, query, &query_layout);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query.
 * @param query_layout The query layout to be retrieved.
 * @return `TILEDB_OK` upon success, and `TILEDB_ERR` upon error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_layout(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    tiledb_layout_t* query_layout) TILEDB_NOEXCEPT;

/**
 * Retrieves the query array.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_query_get_array(ctx, query, &array);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query.
 * @param array The query array to be retrieved.
 * @return `TILEDB_OK` upon success, and `TILEDB_ERR` upon error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_array(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    tiledb_array_t** array) TILEDB_NOEXCEPT;
/**
 * Adds a 1D range along a subarray dimension index, which is in the form
 * (start, end, stride). The datatype of the range components
 * must be the same as the type of the domain of the array in the query.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t dim_idx = 2;
 * int64_t start = 10;
 * int64_t end = 20;
 * tiledb_query_add_range(ctx, query, dim_idx, &start, &end, nullptr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query to add the range to.
 * @param dim_idx The index of the dimension to add the range to.
 * @param start The range start.
 * @param end The range end.
 * @param stride The range stride.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note The stride is currently unsupported. Use `nullptr` as the
 *     stride argument.
 */

TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_add_range(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    uint32_t dim_idx,
    const void* start,
    const void* end,
    const void* stride) TILEDB_NOEXCEPT;

/**
 * Adds a 1D range along a subarray dimension name, which is in the form
 * (start, end, stride). The datatype of the range components
 * must be the same as the type of the domain of the array in the query.
 *
 * **Example:**
 *
 * @code{.c}
 * char* dim_name = "rows";
 * int64_t start = 10;
 * int64_t end = 20;
 * tiledb_query_add_range_by_name(ctx, query, dim_name, &start, &end, nullptr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query to add the range to.
 * @param dim_name The name of the dimension to add the range to.
 * @param start The range start.
 * @param end The range end.
 * @param stride The range stride.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note The stride is currently unsupported. Use `nullptr` as the
 *     stride argument.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_add_range_by_name(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* dim_name,
    const void* start,
    const void* end,
    const void* stride) TILEDB_NOEXCEPT;

/**
 * Adds a 1D variable-sized range along a subarray dimension index, which is in
 * the form (start, end). Applicable only to variable-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t dim_idx = 2;
 * char start[] = "a";
 * char end[] = "bb";
 * tiledb_query_add_range_var(ctx, query, dim_idx, start, 1, end, 2);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query to add the range to.
 * @param dim_idx The index of the dimension to add the range to.
 * @param start The range start.
 * @param start_size The size of the range start in bytes.
 * @param end The range end.
 * @param end_size The size of the range end in bytes.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_add_range_var(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    uint32_t dim_idx,
    const void* start,
    uint64_t start_size,
    const void* end,
    uint64_t end_size) TILEDB_NOEXCEPT;

/**
 * Adds a 1D variable-sized range along a subarray dimension name, which is in
 * the form (start, end). Applicable only to variable-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * char* dim_name = "rows";
 * char start[] = "a";
 * char end[] = "bb";
 * tiledb_query_add_range_var_by_name(ctx, query, dim_name, start, 1, end, 2);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The query to add the range to.
 * @param dim_name The name of the dimension to add the range to.
 * @param start The range start.
 * @param start_size The size of the range start in bytes.
 * @param end The range end.
 * @param end_size The size of the range end in bytes.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_add_range_var_by_name(
    tiledb_ctx_t* ctx,
    tiledb_query_t* query,
    const char* dim_name,
    const void* start,
    uint64_t start_size,
    const void* end,
    uint64_t end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of ranges of the query subarray along a given dimension
 * index.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t range_num;
 * tiledb_query_get_range_num(ctx, query, dim_idx, &range_num);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param dim_idx The index of the dimension whose range number to retrieve.
 * @param range_num The number of ranges to retrieve.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_range_num(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    uint32_t dim_idx,
    uint64_t* range_num) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of ranges of the query subarray along a given dimension
 * name.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t range_num;
 * tiledb_query_get_range_num_from_name(ctx, query, dim_name, &range_num);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param dim_name The name of the dimension whose range number to retrieve.
 * @param range_num The number of ranges to retrieve.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_range_num_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    const char* dim_name,
    uint64_t* range_num) TILEDB_NOEXCEPT;

/**
 * Retrieves a specific range of the query subarray along a given dimension
 * index.
 *
 * **Example:**
 *
 * @code{.c}
 * const void* start;
 * const void* end;
 * const void* stride;
 * tiledb_query_get_range(
 *     ctx, query, dim_idx, range_idx, &start, &end, &stride);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param dim_idx The index of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start The range start to retrieve.
 * @param end The range end to retrieve.
 * @param stride The range stride to retrieve.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_range(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    uint32_t dim_idx,
    uint64_t range_idx,
    const void** start,
    const void** end,
    const void** stride) TILEDB_NOEXCEPT;

/**
 * Retrieves a specific range of the query subarray along a given dimension
 * name.
 *
 * **Example:**
 *
 * @code{.c}
 * const void* start;
 * const void* end;
 * const void* stride;
 * tiledb_query_get_range_from_name(
 *     ctx, query, dim_name, range_idx, &start, &end, &stride);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param dim_name The name of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start The range start to retrieve.
 * @param end The range end to retrieve.
 * @param stride The range stride to retrieve.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_range_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    const char* dim_name,
    uint64_t range_idx,
    const void** start,
    const void** end,
    const void** stride) TILEDB_NOEXCEPT;

/**
 * Retrieves a range's start and end size for a given variable-length
 * dimension index at a given range index.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t start_size;
 * uint64_t end_size;
 * tiledb_query_get_range_var_size(
 *     ctx, query, dim_idx, range_idx, &start_size, &end_size);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param dim_idx The index of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start_size range start size in bytes
 * @param end_size range end size in bytes
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_range_var_size(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    uint32_t dim_idx,
    uint64_t range_idx,
    uint64_t* start_size,
    uint64_t* end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves a range's start and end size for a given variable-length
 * dimension name at a given range index.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t start_size;
 * uint64_t end_size;
 * tiledb_query_get_range_var_size_from_name(
 *     ctx, query, dim_name, range_idx, &start_size, &end_size);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param dim_name The name of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start_size range start size in bytes
 * @param end_size range end size in bytes
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_range_var_size_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    const char* dim_name,
    uint64_t range_idx,
    uint64_t* start_size,
    uint64_t* end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves a specific range of the query subarray along a given
 * variable-length dimension index.
 *
 * **Example:**
 *
 * @code{.c}
 * const void* start;
 * const void* end;
 * tiledb_query_get_range_var(
 *     ctx, query, dim_idx, range_idx, &start, &end);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param dim_idx The index of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start The range start to retrieve.
 * @param end The range end to retrieve.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_range_var(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    uint32_t dim_idx,
    uint64_t range_idx,
    void* start,
    void* end) TILEDB_NOEXCEPT;

/**
 * Retrieves a specific range of the query subarray along a given
 * variable-length dimension name.
 *
 * **Example:**
 *
 * @code{.c}
 * const void* start;
 * const void* end;
 * tiledb_query_get_range_var_from_name(
 *     ctx, query, dim_name, range_idx, &start, &end);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param dim_name The name of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start The range start to retrieve.
 * @param end The range end to retrieve.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_query_get_range_var_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    const char* dim_name,
    uint64_t range_idx,
    void* start,
    void* end) TILEDB_NOEXCEPT;

/**
 * Retrieves the estimated result size for a fixed-sized attribute/dimension.
 * This is an estimate and may not be sufficient to read all results for the
 * requested range, in particular for sparse arrays or array with
 * var-length attributes.
 * Query status must be checked and resubmitted if not complete.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t size;
 * tiledb_query_get_est_result_size(ctx, query, "a", &size);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param name The attribute/dimension name.
 * @param size The size (in bytes) to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_est_result_size(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    const char* name,
    uint64_t* size) TILEDB_NOEXCEPT;

/**
 * Retrieves the estimated result size for a var-sized attribute/dimension.
 * This is an estimate and may not be sufficient to read all results for the
 * requested range, for sparse arrays or any array with
 * var-length attributes.
 * Query status must be checked and resubmitted if not complete.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t size_off, size_val;
 * tiledb_query_get_est_result_size_var(
 *     ctx, query, "a", &size_off, &size_val);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param name The attribute/dimension name.
 * @param size_off The size of the offsets (in bytes) to be retrieved.
 * @param size_val The size of the values (in bytes) to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_est_result_size_var(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    const char* name,
    uint64_t* size_off,
    uint64_t* size_val) TILEDB_NOEXCEPT;

/**
 * Retrieves the estimated result size for a fixed-sized, nullable attribute.
 * This is an estimate and may not be sufficient to read all results for the
 * requested range, for sparse arrays or any array with
 * var-length attributes.
 * Query status must be checked and resubmitted if not complete.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t size_val;
 * uint64_t size_validity;
 * tiledb_query_get_est_result_size_nullable(ctx, query, "a", &size_val,
 * &size_validity);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param name The attribute name.
 * @param size_val The size of the values (in bytes) to be retrieved.
 * @param size_validity The size of the validity values (in bytes) to be
 * retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_est_result_size_nullable(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    const char* name,
    uint64_t* size_val,
    uint64_t* size_validity) TILEDB_NOEXCEPT;

/**
 * Retrieves the estimated result size for a var-sized, nullable attribute.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t size_off, size_val, size_validity;
 * tiledb_query_get_est_result_size_var_nullable(
 *     ctx, query, "a", &size_off, &size_val, &size_validity);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param name The attribute name.
 * @param size_off The size of the offsets (in bytes) to be retrieved.
 * @param size_val The size of the values (in bytes) to be retrieved.
 * @param size_validity The size of the validity values (in bytes) to be
 * retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_est_result_size_var_nullable(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    const char* name,
    uint64_t* size_off,
    uint64_t* size_val,
    uint64_t* size_validity) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of written fragments. Applicable only to WRITE
 * queries.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t num;
 * tiledb_query_get_fragment_num(ctx, query, &num);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param num The number of written fragments to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_fragment_num(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    uint32_t* num) TILEDB_NOEXCEPT;

/**
 * Retrieves the URI of the written fragment with the input index. Applicable
 * only to WRITE queries.
 *
 * **Example:**
 *
 * @code{.c}
 * const char* uri;
 * tiledb_query_get_fragment_uri(
 *     ctx, query, 0, &uri);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param idx The index of the written fragment.
 * @param uri The URI of the written fragment to be returned.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note Make sure to make a copy of `uri` after its retrieval, as the
 *     constant pointer may be updated internally as new fragments
 *     are being written.
 */
TILEDB_EXPORT int32_t tiledb_query_get_fragment_uri(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    uint64_t idx,
    const char** uri) TILEDB_NOEXCEPT;

/**
 * Retrieves the timestamp range of the written fragment with the input index.
 * Applicable only to WRITE queries.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t t1, t2;
 * tiledb_query_get_fragment_timestamp_range(
 *     ctx, query, 0, &t1, &t2);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param query The query.
 * @param idx The index of the written fragment.
 * @param t1 The start value of the timestamp range of the
 *     written fragment to be returned.
 * @param t2 The end value of the timestamp range of the
 *     written fragment to be returned.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_get_fragment_timestamp_range(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    uint64_t idx,
    uint64_t* t1,
    uint64_t* t2) TILEDB_NOEXCEPT;

/**
 * Return a TileDB subarray object from the given query.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_subarray_t* subarray;
 * tiledb_query_get_subarray_t(array, &subarray);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array An open array object.
 * @param subarray The retrieved subarray object if available.
 * @return `TILEDB_OK` for success or `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT
int32_t tiledb_query_get_subarray_t(
    tiledb_ctx_t* ctx,
    const tiledb_query_t* query,
    tiledb_subarray_t** subarray) TILEDB_NOEXCEPT;

/* ****************************** */
/*          QUERY CONDITION       */
/* ****************************** */

/**
 * Allocates a TileDB query condition object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_condition_t* query_condition;
 * tiledb_query_condition_alloc(ctx, &query_condition);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param cond The allocated query condition object.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_condition_alloc(
    tiledb_ctx_t* ctx, tiledb_query_condition_t** cond) TILEDB_NOEXCEPT;

/**
 * Frees a TileDB query condition object.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t value = 5;
 * tiledb_query_condition_t* query_condition;
 * tiledb_query_condition_alloc(
 *   ctx, "longitude", &value, sizeof(value), TILEDB_LT, &query_condition);
 * tiledb_query_set_condition(ctx, query, query_condition);
 * tiledb_query_submit(ctx, query);
 * tiledb_query_condition_free(&query_condition);
 * @endcode
 *
 * @param cond The query condition object to be freed.
 */
TILEDB_EXPORT void tiledb_query_condition_free(tiledb_query_condition_t** cond)
    TILEDB_NOEXCEPT;

/**
 * Initializes a TileDB query condition object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_condition_t* query_condition;
 * tiledb_query_condition_alloc(ctx, &query_condition);
 *
 * uint32_t value = 5;
 * tiledb_query_condition_init(
 *   ctx, query_condition, "longitude", &value, sizeof(value), TILEDB_LT);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param cond The allocated query condition object.
 * @param attribute_name The attribute name.
 * @param condition_value The value to compare against an attribute value.
 * @param condition_value_size The byte size of `condition_value`.
 * @param op The comparison operator.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_condition_init(
    tiledb_ctx_t* ctx,
    tiledb_query_condition_t* cond,
    const char* attribute_name,
    const void* condition_value,
    uint64_t condition_value_size,
    tiledb_query_condition_op_t op) TILEDB_NOEXCEPT;

/**
 * Combines two query condition objects into a newly allocated
 * condition. Does not mutate or free the input condition objects.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_query_condition_t* query_condition_1;
 * tiledb_query_condition_alloc(ctx, &query_condition_1);
 * uint32_t value_1 = 5;
 * tiledb_query_condition_init(
 *   ctx,
 *   query_condition_1,
 *   "longitude",
 *   &value_1,
 *   sizeof(value_1),
 *   TILEDB_LT);
 *
 * tiledb_query_condition_t* query_condition_2;
 * tiledb_query_condition_alloc(ctx, &query_condition_2);
 * uint32_t value_2 = 20;
 * tiledb_query_condition_init(
 *   ctx,
 *   query_condition_2,
 *   "latitude",
 *   &value_2,
 *   sizeof(value_2),
 *   TILEDB_GE);
 *
 * tiledb_query_condition_t* query_condition_3;
 * tiledb_query_condition_combine(
 *   ctx, query_condition_1, query_condition_2, TILEDB_AND, &query_condition_3);
 *
 * tiledb_query_condition_free(&query_condition_1);
 * tiledb_query_condition_free(&query_condition_2);
 *
 * tiledb_query_set_condition(ctx, query, query_condition_3);
 * tiledb_query_submit(ctx, query);
 * tiledb_query_condition_free(&query_condition_3);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param attribute_name The attribute name.
 * @param condition_value The value to compare against an attribute value.
 * @param condition_value_size The byte size of `condition_value`.
 * @param op The comparison operator.
 * @param cond The allocated query condition object.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_query_condition_combine(
    tiledb_ctx_t* ctx,
    const tiledb_query_condition_t* left_cond,
    const tiledb_query_condition_t* right_cond,
    tiledb_query_condition_combination_op_t combination_op,
    tiledb_query_condition_t** combined_cond) TILEDB_NOEXCEPT;

/* ********************************* */
/*             SUBARRAY              */
/* ********************************* */

/**
 * Allocates a TileDB subarray object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_subarray_t* subarray;
 * tiledb_subarray_alloc(ctx, array, &subarray);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array An open array object.
 * @param subarray The subarray object to be created.
 * @return `TILEDB_OK` for success or `TILEDB_OOM` or `TILEDB_ERR` for error.
 *
 * @note The allocated subarray initially has internal coalesce_ranges == true.
 */
TILEDB_EXPORT int32_t tiledb_subarray_alloc(
    tiledb_ctx_t* ctx,
    const tiledb_array_t* array,
    tiledb_subarray_t** subarray) TILEDB_NOEXCEPT;

/**
 * Set the subarray config.
 *
 * Setting the configuration with this function overrides the following
 * Subarray-level parameters only:
 *
 * - `sm.read_range_oob`
 */
TILEDB_EXPORT int32_t tiledb_subarray_set_config(
    tiledb_ctx_t* ctx,
    tiledb_subarray_t* subarray,
    tiledb_config_t* config) TILEDB_NOEXCEPT;

/**
 * Frees a TileDB subarray object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_subarray_t* subarray;
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * tiledb_subarray_alloc(ctx, array, &subarray);
 * tiledb_array_close(ctx, array);
 * tiledb_subarray_free(&subarray);
 * @endcode
 *
 * @param subarray The subarray object to be freed.
 */
TILEDB_EXPORT void tiledb_subarray_free(tiledb_subarray_t** subarray)
    TILEDB_NOEXCEPT;

/**
 * Set coalesce_ranges property on a TileDB subarray object.
 * Intended to be used just after tiledb_subarray_alloc() to replace
 * the initial coalesce_ranges == true
 * with coalesce_ranges = false if
 * needed.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_subarray_t* subarray;
 * //tiledb_subarray_alloc internally defaults to 'coalesce_ranges == true'
 * tiledb_subarray_alloc(ctx, array, &subarray);
 * // so manually set to 'false' to match earlier behaviour with older
 * // tiledb_query_ subarray actions.
 * bool coalesce_ranges = false;
 * tiledb_subarray_set_coalesce_ranges(ctx, subarray, coalesce_ranges);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param subarray The subarray object to change.
 * @param coalesce_ranges The true/false value to be set
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_set_coalesce_ranges(
    tiledb_ctx_t* ctx,
    tiledb_subarray_t* subarray,
    int coalesce_ranges) TILEDB_NOEXCEPT;

/**
 * Populates a subarray with specific indicies.
 *
 * **Example:**
 *
 * The following sets a 2D subarray [0,10], [20, 30] to the subarray.
 *
 * @code{.c}
 * tiledb_subarray_t *subarray;
 * uint64_t subarray_v[] = { 0, 10, 20, 30};
 * tiledb_subarray_set_subarray(ctx, subarray, subarray_v);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param subarray The TileDB subarray object.
 * @param subarray_v The subarray values which can be used to limit the subarray
 * read/write.
 *     It should be a sequence of [low, high] pairs (one pair per dimension).
 *     When the subarray is used for writes, this is meaningful only
 *     for dense arrays, and specifically dense writes. Note that `subarray_a`
 *     must have the same type as the domain of the subarray's associated
 *     array.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT
int32_t tiledb_subarray_set_subarray(
    tiledb_ctx_t* ctx,
    tiledb_subarray_t* subarray_s,
    const void* subarray_v) TILEDB_NOEXCEPT;

/**
 * Adds a 1D range along a subarray dimension index, which is in the form
 * (start, end, stride). The datatype of the range components
 * must be the same as the type of the domain of the array in the query.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t dim_idx = 2;
 * int64_t start = 10;
 * int64_t end = 20;
 * tiledb_subarray_add_range(ctx, subarray, dim_idx, &start, &end, nullptr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param query The subarray to add the range to.
 * @param dim_idx The index of the dimension to add the range to.
 * @param start The range start.
 * @param end The range end.
 * @param stride The range stride.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 *
 * @note The stride is currently unsupported. Use 0/NULL/nullptr as the
 *     stride argument.
 */
TILEDB_EXPORT int32_t tiledb_subarray_add_range(
    tiledb_ctx_t* ctx,
    tiledb_subarray_t* subarray,
    uint32_t dim_idx,
    const void* start,
    const void* end,
    const void* stride) TILEDB_NOEXCEPT;

/**
 * Adds a 1D range along a subarray dimension name, which is in the form
 * (start, end, stride). The datatype of the range components
 * must be the same as the type of the domain of the array in the query.
 *
 * **Example:**
 *
 * @code{.c}
 * char* dim_name = "rows";
 * int64_t start = 10;
 * int64_t end = 20;
 * tiledb_subarray_add_range_by_name(
 *     ctx, subarray, dim_name, &start, &end, nullptr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param subarray The subarray to add the range to.
 * @param dim_name The name of the dimension to add the range to.
 * @param start The range start.
 * @param end The range end.
 * @param stride The range stride.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 *
 * @note The stride is currently unsupported. Use 0/NULL/nullptr as the
 *     stride argument.
 */
TILEDB_EXPORT int32_t tiledb_subarray_add_range_by_name(
    tiledb_ctx_t* ctx,
    tiledb_subarray_t* subarray,
    const char* dim_name,
    const void* start,
    const void* end,
    const void* stride) TILEDB_NOEXCEPT;

/**
 * Adds a 1D variable-sized range along a subarray dimension index, which is in
 * the form (start, end). Applicable only to variable-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t dim_idx = 2;
 * char start[] = "a";
 * char end[] = "bb";
 * tiledb_subarray_add_range_var(ctx, subarray, dim_idx, start, 1, end, 2);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param subarray The subarray to add the range to.
 * @param dim_idx The index of the dimension to add the range to.
 * @param start The range start.
 * @param start_size The size of the range start in bytes.
 * @param end The range end.
 * @param end_size The size of the range end in bytes.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_add_range_var(
    tiledb_ctx_t* ctx,
    tiledb_subarray_t* subarray,
    uint32_t dim_idx,
    const void* start,
    uint64_t start_size,
    const void* end,
    uint64_t end_size) TILEDB_NOEXCEPT;

/**
 * Adds a 1D variable-sized range along a subarray dimension name, which is in
 * the form (start, end). Applicable only to variable-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * char* dim_name = "rows";
 * char start[] = "a";
 * char end[] = "bb";
 * tiledb_subarray_add_range_var_by_name(
 *     ctx, subarray, dim_name, start, 1, end, 2);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param subarray The subarray to add the range to.
 * @param dim_name The name of the dimension to add the range to.
 * @param start The range start.
 * @param start_size The size of the range start in bytes.
 * @param end The range end.
 * @param end_size The size of the range end in bytes.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_add_range_var_by_name(
    tiledb_ctx_t* ctx,
    tiledb_subarray_t* subarray,
    const char* dim_name,
    const void* start,
    uint64_t start_size,
    const void* end,
    uint64_t end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of ranges of the query subarray along a given dimension
 * index.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t range_num;
 * tiledb_subarray_get_range_num(ctx, subarray, dim_idx, &range_num);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param subarray The subarray.
 * @param dim_idx The index of the dimension for which to retrieve number of
 * ranges.
 * @param range_num Receives the retrieved number of ranges.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_get_range_num(
    tiledb_ctx_t* ctx,
    const tiledb_subarray_t* subarray,
    uint32_t dim_idx,
    uint64_t* range_num) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of ranges of the subarray along a given dimension
 * name.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t range_num;
 * tiledb_subarray_get_range_num_from_name(ctx, subarray, dim_name, &range_num);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param subarray The subarray.
 * @param dim_name The name of the dimension whose range number to retrieve.
 * @param range_num Receives the retrieved number of ranges.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_get_range_num_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_subarray_t* subarray,
    const char* dim_name,
    uint64_t* range_num) TILEDB_NOEXCEPT;

/**
 * Retrieves a specific range of the subarray along a given dimension
 * index.
 *
 * **Example:**
 *
 * @code{.c}
 * const void* start;
 * const void* end;
 * const void* stride;
 * tiledb_subarray_get_range(
 *     ctx, subarray, dim_idx, range_idx, &start, &end, &stride);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param subarray The subarray.
 * @param dim_idx The index of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start Receives the retrieved range start.
 * @param end Receives the received range end.
 * @param stride Receives the retrieved range stride.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_get_range(
    tiledb_ctx_t* ctx,
    const tiledb_subarray_t* subarray,
    uint32_t dim_idx,
    uint64_t range_idx,
    const void** start,
    const void** end,
    const void** stride) TILEDB_NOEXCEPT;

/**
 * Retrieves a specific range of the subarray along a given dimension
 * name.
 *
 * **Example:**
 *
 * @code{.c}
 * const void* start;
 * const void* end;
 * const void* stride;
 * tiledb_subarray_get_range_from_name(
 *     ctx, query, dim_name, range_idx, &start, &end, &stride);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param subarray The subarray.
 * @param dim_name The name of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start Receives the retrieved range start.
 * @param end Receives the retrieved range end.
 * @param stride Receives the retrieved range stride.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_get_range_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_subarray_t* subarray,
    const char* dim_name,
    uint64_t range_idx,
    const void** start,
    const void** end,
    const void** stride) TILEDB_NOEXCEPT;

/**
 * Retrieves a range's start and end size for a given variable-length
 * dimension index at a given range index.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t start_size;
 * uint64_t end_size;
 * tiledb_subarray_get_range_var_size(
 *     ctx, subarray, dim_idx, range_idx, &start_size, &end_size);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param subarray The subarray.
 * @param dim_idx The index of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start_size Receives the retrieved range start size in bytes
 * @param end_size Receives the retrieved range end size in bytes
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_get_range_var_size(
    tiledb_ctx_t* ctx,
    const tiledb_subarray_t* subarray,
    uint32_t dim_idx,
    uint64_t range_idx,
    uint64_t* start_size,
    uint64_t* end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves a range's start and end size for a given variable-length
 * dimension name at a given range index.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t start_size;
 * uint64_t end_size;
 * tiledb_subarray_get_range_var_size_from_name(
 *     ctx, subarray, dim_name, range_idx, &start_size, &end_size);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param subarray The subarray.
 * @param dim_name The name of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start_size Receives the retrieved range start size in bytes
 * @param end_size Receives the retrieved range end size in bytes
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_get_range_var_size_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_subarray_t* subarray,
    const char* dim_name,
    uint64_t range_idx,
    uint64_t* start_size,
    uint64_t* end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves a specific range of the subarray along a given
 * variable-length dimension index.
 *
 * **Example:**
 *
 * @code{.c}
 * const void* start;
 * const void* end;
 * tiledb_subarray_get_range_var(
 *     ctx, subarray, dim_idx, range_idx, &start, &end);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param subarray The subarray.
 * @param dim_idx The index of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start Receives the retrieved range start.
 * @param end Receives the retrieved range end.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_get_range_var(
    tiledb_ctx_t* ctx,
    const tiledb_subarray_t* subarray,
    uint32_t dim_idx,
    uint64_t range_idx,
    void* start,
    void* end) TILEDB_NOEXCEPT;

/**
 * Retrieves a specific range of the subarray along a given
 * variable-length dimension name.
 *
 * **Example:**
 *
 * @code{.c}
 * const void* start;
 * const void* end;
 * tiledb_subarray_get_range_var_from_name(
 *     ctx, subarray, dim_name, range_idx, &start, &end);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param subarray The subarray.
 * @param dim_name The name of the dimension to retrieve the range from.
 * @param range_idx The index of the range to retrieve.
 * @param start Receives the retrieved range start.
 * @param end Receives the retrieved range end.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_subarray_get_range_var_from_name(
    tiledb_ctx_t* ctx,
    const tiledb_subarray_t* subarray,
    const char* dim_name,
    uint64_t range_idx,
    void* start,
    void* end) TILEDB_NOEXCEPT;

/* ********************************* */
/*               ARRAY               */
/* ********************************* */

/**
 * Allocates a TileDB array object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The array URI.
 * @param array The array object to be created.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_alloc(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_array_t** array) TILEDB_NOEXCEPT;

/**
 * Sets the starting timestamp to use when opening (and reopening) the array.
 * This is an inclusive bound. The default value is `0`.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "s3://tiledb_bucket/my_array", &array);
 * tiledb_array_set_open_timestamp_start(ctx, array, 1234);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array to set the timestamp on.
 * @param timestamp_start The epoch timestamp in milliseconds.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_set_open_timestamp_start(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint64_t timestamp_start) TILEDB_NOEXCEPT;

/**
 * Sets the ending timestamp to use when opening (and reopening) the array.
 * This is an inclusive bound. The UINT64_MAX timestamp is a reserved timestamp
 * that will be interpretted as the current timestamp when an array is opened.
 * The default value is `UINT64_MAX`.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "s3://tiledb_bucket/my_array", &array);
 * tiledb_array_set_open_timestamp_end(ctx, array, 5678);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array to set the timestamp on.
 * @param timestamp_end The epoch timestamp in milliseconds. Use UINT64_MAX for
 *   the current timestamp.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_set_open_timestamp_end(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint64_t timestamp_end) TILEDB_NOEXCEPT;

/**
 * Gets the starting timestamp used when opening (and reopening) the array.
 * This is an inclusive bound.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "s3://tiledb_bucket/my_array", &array);
 * tiledb_array_set_open_timestamp_start(ctx, array, 1234);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 *
 * uint64_t timestamp_start;
 * tiledb_array_get_open_timestamp_start(ctx, array, &timestamp_start);
 * assert(timestamp_start == 1234);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array to set the timestamp on.
 * @param timestamp_start The output epoch timestamp in milliseconds.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_open_timestamp_start(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint64_t* timestamp_start) TILEDB_NOEXCEPT;

/**
 * Gets the ending timestamp used when opening (and reopening) the array.
 * This is an inclusive bound. If UINT64_MAX was set, this will return
 * the timestamp at the time the array was opened. If the array has not
 * yet been opened, it will return UINT64_MAX.`
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "s3://tiledb_bucket/my_array", &array);
 * tiledb_array_set_open_timestamp_end(ctx, array, 5678);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 *
 * uint64_t timestamp_end;
 * tiledb_array_get_open_timestamp_end(ctx, array, &timestamp_end);
 * assert(timestamp_start == 5678);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array to set the timestamp on.
 * @param timestamp_end The output epoch timestamp in milliseconds.
 * @return `TILEDB_OK` for success or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_open_timestamp_end(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint64_t* timestamp_end) TILEDB_NOEXCEPT;

/**
 * Deletes all written array data.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_delete_array(
 *   ctx, array, "hdfs:///temp/my_array");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array to delete the data from.
 * @param uri The Array's URI.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_delete_array(
    tiledb_ctx_t* ctx, tiledb_array_t* array, const char* uri) TILEDB_NOEXCEPT;

/**
 * Deletes array fragments written between the input timestamps.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_delete_fragments(
 *   ctx, array, "hdfs:///temp/my_array", 0, UINT64_MAX);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array to delete the fragments from.
 * @param uri The URI of the fragments' parent Array.
 * @param timestamp_start The epoch timestamp in milliseconds.
 * @param timestamp_end The epoch timestamp in milliseconds. Use UINT64_MAX for
 *   the current timestamp.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_delete_fragments(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    const char* uri,
    uint64_t timestamp_start,
    uint64_t timestamp_end) TILEDB_NOEXCEPT;

/**
 * Opens a TileDB array. The array is opened using a query type as input.
 * This is to indicate that queries created for this `tiledb_array_t`
 * object will inherit the query type. In other words, `tiledb_array_t`
 * objects are opened to receive only one type of queries.
 * They can always be closed and be re-opened with another query type.
 * Also there may be many different `tiledb_array_t`
 * objects created and opened with different query types.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array object to be opened.
 * @param query_type The type of queries the array object will be receiving.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note If the same array object is opened again without being closed,
 *     an error will be thrown.
 * @note The config should be set before opening an array.
 * @note If the array is to be opened at a specfic time interval, the
 *      `timestamp{start, end}` values should be set to a config that's set to
 *       the array object before opening the array.
 */
TILEDB_EXPORT int32_t tiledb_array_open(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    tiledb_query_type_t query_type) TILEDB_NOEXCEPT;

/**
 * Similar to `tiledb_array_open`, but this function takes as input a
 * timestamp, representing time in milliseconds ellapsed since
 * 1970-01-01 00:00:00 +0000 (UTC). Opening the array at a
 * timestamp provides a view of the array with all writes/updates that
 * happened at or before `timestamp` (i.e., excluding those that
 * occurred after `timestamp`). This function is useful to ensure
 * consistency at a potential distributed setting, where machines
 * need to operate on the same view of the array.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * // Assuming `timestamp` is time represented in milliseconds:
 * tiledb_array_open_at(ctx, array, TILEDB_READ, timestamp);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array object to be opened.
 * @param query_type The type of queries the array object will be receiving.
 * @param timestamp The timestamp to open the array at.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note If the same array object is opened again without being closed,
 *     an error will be thrown.
 * @note This function is applicable only to read queries.
 * @note The config should be set before opening an array.
 * @note If the array is to be opened at a specfic time interval, the
 *      `timestamp{start, end}` values should be set to a config that's set to
 *       the array object before opening the array.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_open_at(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    tiledb_query_type_t query_type,
    uint64_t timestamp) TILEDB_NOEXCEPT;

/**
 * Opens an encrypted array using the given encryption key. This function has
 * the same semantics as `tiledb_array_open()` but is used for encrypted arrays.
 *
 * An encrypted array must be opened with this function before queries can be
 * issued to it.
 *
 * **Example:**
 *
 * @code{.c}
 * // Load AES-256 key from disk, environment variable, etc.
 * uint8_t key[32] = ...;
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * tiledb_array_open_with_key(ctx, array, TILEDB_READ,
 *     TILEDB_AES_256_GCM, key, sizeof(key));
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array object to be opened.
 * @param query_type The type of queries the array object will be receiving.
 * @param encryption_type The encryption type to use.
 * @param encryption_key The encryption key to use.
 * @param key_length Length in bytes of the encryption key.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note The config should be set before opening an array.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_open_with_key(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    tiledb_query_type_t query_type,
    tiledb_encryption_type_t encryption_type,
    const void* encryption_key,
    uint32_t key_length) TILEDB_NOEXCEPT;

/**
 * Similar to `tiledb_array_open_with_key`, but this function takes as
 * input a timestamp, representing time in milliseconds ellapsed since
 * 1970-01-01 00:00:00 +0000 (UTC). Opening the array at a
 * timestamp provides a view of the array with all writes/updates that
 * happened at or before `timestamp` (i.e., excluding those that
 * occurred after `timestamp`). This function is useful to ensure
 * consistency at a potential distributed setting, where machines
 * need to operate on the same view of the array.
 *
 * **Example:**
 *
 * @code{.c}
 * // Load AES-256 key from disk, environment variable, etc.
 * uint8_t key[32] = ...;
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * // Assuming `timestamp` is time represented in milliseconds:
 * tiledb_array_open_at_with_key(ctx, array, TILEDB_READ,
 *     TILEDB_AES_256_GCM, key, sizeof(key), timestamp);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array object to be opened.
 * @param query_type The type of queries the array object will be receiving.
 * @param encryption_type The encryption type to use.
 * @param encryption_key The encryption key to use.
 * @param key_length Length in bytes of the encryption key.
 * @param timestamp The timestamp to open the array at.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note If the same array object is opened again without being closed,
 *     an error will be thrown.
 * @note This function is applicable only to read queries.
 * @note The config should be set before opening an array.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_open_at_with_key(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    tiledb_query_type_t query_type,
    tiledb_encryption_type_t encryption_type,
    const void* encryption_key,
    uint32_t key_length,
    uint64_t timestamp) TILEDB_NOEXCEPT;

/**
 * Checks if the array is open.
 *
 * @param ctx The TileDB context.
 * @param array The array to be checked.
 * @param is_open `1` if the array is open and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_is_open(
    tiledb_ctx_t* ctx, tiledb_array_t* array, int32_t* is_open) TILEDB_NOEXCEPT;

/**
 * Reopens a TileDB array (the array must be already open). This is useful
 * when the array got updated after it got opened and the `tiledb_array_t`
 * object got created. To sync-up with the updates, the user must either
 * close the array and open with `tiledb_array_open`, or just use
 * `tiledb_array_reopen` without closing. This function will be generally
 * faster than the former alternative.
 *
 * Note: reopening encrypted arrays does not require the encryption key.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * tiledb_array_reopen(ctx, array);
 *
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array object to be re-opened.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note This is applicable only to arrays opened for reads.
 * @note If the array is to be reopened after opening at a specfic time
 *      interval, the `timestamp{start, end}` values and subsequent config
 *      object should be reset for the array before reopening.
 */
TILEDB_EXPORT int32_t
tiledb_array_reopen(tiledb_ctx_t* ctx, tiledb_array_t* array) TILEDB_NOEXCEPT;

/**
 * Reopens a TileDB array (the array must be already open) at a specific
 * timestamp.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * uint64_t timestamp = tiledb_timestamp_now_ms();
 * tiledb_array_reopen_at(ctx, array, timestamp);
 *
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array object to be re-opened.
 * @param timestamp Timestamp at which to reopen.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note This is applicable only to arrays opened for reads.
 * @note If the array is to be reopened after opening at a specfic time
 *      interval, the `timestamp{start, end}` values and subsequent config
 *      object should be reset for the array before reopening.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_reopen_at(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint64_t timestamp) TILEDB_NOEXCEPT;

/**
 * The start/end timestamps for opening an array
 * are now set in the config.
 *
 * Returns the timestamp, representing time in milliseconds ellapsed since
 * 1970-01-01 00:00:00 +0000 (UTC), at which the array was opened. See also the
 * documentation of `tiledb_array_open_at`.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "s3://tiledb_bucket/my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * // Get the timestamp the array at which the array was opened.
 * uint64_t timestamp;
 * tiledb_array_get_timestamp(ctx, array, &timestamp);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array to retrieve the timestamp for.
 * @param timestamp Set to the timestamp at which the array was opened.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note The array does not need to be open to use this function.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_get_timestamp(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint64_t* timestamp) TILEDB_NOEXCEPT;

/**
 * Sets the array config.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "s3://tiledb_bucket/my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * // Set the config for the given array.
 * tiledb_config_t* config;
 * tiledb_array_set_config(ctx, array, config);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array to set the config for.
 * @param config The config to be set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note The array does not need to be opened via `tiledb_array_open_at` to use
 *      this function.
 * @note The config should be set before opening an array.
 */
TILEDB_EXPORT int32_t tiledb_array_set_config(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    tiledb_config_t* config) TILEDB_NOEXCEPT;

/**
 * Gets the array config.
 *
 * **Example:**
 *
 * @code{.c}
 * // Retrieve the config for the given array.
 * tiledb_config_t* config;
 * tiledb_array_get_config(ctx, array, config);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array to set the config for.
 * @param config Set to the retrieved config.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_config(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    tiledb_config_t** config) TILEDB_NOEXCEPT;

/**
 * Closes a TileDB array.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * tiledb_array_close(ctx, array);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array object to be closed.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note If the array object has already been closed, the function has
 *     no effect.
 */
TILEDB_EXPORT int32_t
tiledb_array_close(tiledb_ctx_t* ctx, tiledb_array_t* array) TILEDB_NOEXCEPT;

/**
 * Frees a TileDB array object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * tiledb_array_close(ctx, array);
 * tiledb_array_free(&array);
 * @endcode
 *
 * @param array The array object to be freed.
 */
TILEDB_EXPORT void tiledb_array_free(tiledb_array_t** array) TILEDB_NOEXCEPT;

/**
 * Retrieves the schema of an array.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_t* array_schema;
 * tiledb_array_get_schema(ctx, array, &array_schema);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The open array.
 * @param array_schema The array schema to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 *
 * @note The user must free the array schema with `tiledb_array_schema_free`.
 */
TILEDB_EXPORT int32_t tiledb_array_get_schema(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    tiledb_array_schema_t** array_schema) TILEDB_NOEXCEPT;

/**
 * Retrieves the query type with which the array was opened.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "hdfs:///tiledb_arrays/my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * tiledb_query_type_t query_type;
 * tiledb_array_get_type(ctx, array, &query_type);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array The array.
 * @param query_type The query type to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_query_type(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    tiledb_query_type_t* query_type) TILEDB_NOEXCEPT;

/**
 * Creates a new TileDB array given an input schema.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_create(ctx, "hdfs:///tiledb_arrays/my_array", array_schema);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The array name.
 * @param array_schema The array schema.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_create(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    const tiledb_array_schema_t* array_schema) TILEDB_NOEXCEPT;

/**
 * Creates a new encrypted TileDB array given an input schema.
 *
 * Encrypted arrays can only be created through this function.
 *
 * **Example:**
 *
 * @code{.c}
 * uint8_t key[32] = ...;
 * tiledb_array_create_with_key(
 *     ctx, "hdfs:///tiledb_arrays/my_array", array_schema,
 *     TILEDB_AES_256_GCM, key, sizeof(key));
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The array name.
 * @param array_schema The array schema.
 * @param encryption_type The encryption type to use.
 * @param encryption_key The encryption key to use.
 * @param key_length Length in bytes of the encryption key.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_create_with_key(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    const tiledb_array_schema_t* array_schema,
    tiledb_encryption_type_t encryption_type,
    const void* encryption_key,
    uint32_t key_length) TILEDB_NOEXCEPT;

/**
 * Depending on the consoliation mode in the config, consolidates either the
 * fragment files, fragment metadata files, or array metadata files into a
 * single file.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_consolidate(
 *     ctx, "hdfs:///tiledb_arrays/my_array", nullptr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The name of the TileDB array whose metadata will
 *     be consolidated.
 * @param config Configuration parameters for the consolidation
 *     (`nullptr` means default, which will use the config from `ctx`).
 *     The `sm.consolidation.mode` parameter determines which type of
 *     consolidation to perform.
 *
 * @return `TILEDB_OK` on success, and `TILEDB_ERR` on error.
 */
TILEDB_EXPORT int32_t tiledb_array_consolidate(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_config_t* config) TILEDB_NOEXCEPT;

/**
 * Depending on the consoliation mode in the config, consolidates either the
 * fragment files, fragment metadata files, or array metadata files into a
 * single file.
 *
 * **Example:**
 *
 * @code{.c}
 * uint8_t key[32] = ...;
 * tiledb_array_consolidate_with_key(
 *     ctx, "hdfs:///tiledb_arrays/my_array",
 *     TILEDB_AES_256_GCM, key, sizeof(key), nullptr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The name of the TileDB array to be consolidated.
 * @param encryption_type The encryption type to use.
 * @param encryption_key The encryption key to use.
 * @param key_length Length in bytes of the encryption key.
 * @param config Configuration parameters for the consolidation
 *     (`nullptr` means default, which will use the config from `ctx`).
 *     The `sm.consolidation.mode` parameter determines which type of
 *     consolidation to perform.
 *
 * @return `TILEDB_OK` on success, and `TILEDB_ERR` on error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_consolidate_with_key(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_encryption_type_t encryption_type,
    const void* encryption_key,
    uint32_t key_length,
    tiledb_config_t* config) TILEDB_NOEXCEPT;

/**
 * Cleans up the array, such as consolidated fragments and array metadata.
 * Note that this will coarsen the granularity of time traveling (see docs
 * for more information).
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_vacuum(ctx, "hdfs:///tiledb_arrays/my_array");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The name of the TileDB array to vacuum.
 * @param config Configuration parameters for the vacuuming
 *     (`nullptr` means default, which will use the config from `ctx`).
 * @return `TILEDB_OK` on success, and `TILEDB_ERR` on error.
 */
TILEDB_EXPORT int32_t tiledb_array_vacuum(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_config_t* config) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain from an array. This is the union of the
 * non-empty domains of the array fragments.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t domain[4]; // Assuming a 2D array, 2 [low, high] pairs
 * int32_t is_empty;
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * tiledb_array_get_non_empty_domain(ctx, array, domain, &is_empty);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param array The array object (must be opened beforehand).
 * @param domain The domain to be retrieved.
 * @param is_empty The function sets it to `1` if the non-empty domain is
 *     empty (i.e., the array does not contain any data yet), and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_non_empty_domain(
    tiledb_ctx_t* ctx, tiledb_array_t* array, void* domain, int32_t* is_empty)
    TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain from an array for a given dimension index.
 * This is the union of the non-empty domains of the array fragments on
 * the given dimension.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t domain[2];
 * int32_t is_empty;
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * tiledb_array_get_non_empty_domain_from_index(ctx, array, 0, domain,
 * &is_empty);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param array The array object (must be opened beforehand).
 * @param idx The dimension index, following the order as it was defined
 *      in the domain of the array schema.
 * @param domain The domain to be retrieved.
 * @param is_empty The function sets it to `1` if the non-empty domain is
 *     empty (i.e., the array does not contain any data yet), and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_non_empty_domain_from_index(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint32_t idx,
    void* domain,
    int32_t* is_empty) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain from an array for a given dimension name.
 * This is the union of the non-empty domains of the array fragments on
 * the given dimension.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t domain[2];
 * int32_t is_empty;
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * tiledb_array_get_non_empty_domain_from_name(ctx, array, "d1", domain,
 * &is_empty);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param array The array object (must be opened beforehand).
 * @param name The dimension name.
 * @param domain The domain to be retrieved.
 * @param is_empty The function sets it to `1` if the non-empty domain is
 *     empty (i.e., the array does not contain any data yet), and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_non_empty_domain_from_name(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    const char* name,
    void* domain,
    int32_t* is_empty) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain range sizes from an array for a given
 * dimension index. This is the union of the non-empty domains of the array
 * fragments on the given dimension. Applicable only to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t is_empty;
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * uint64_t start_size, end_size;
 * tiledb_array_get_non_empty_domain_var_size_from_index(
 *     ctx, array, 0, &start_size, &end_size, &is_empty);
 * // If non-empty domain range is `[aa, dddd]`, then `start_size = 2`
 * // and `end_size = 4`.
 * @endcode
 *
 * @param ctx The TileDB context
 * @param array The array object (must be opened beforehand).
 * @param idx The dimension index, following the order as it was defined
 *      in the domain of the array schema.
 * @param start_size The size in bytes of the start range.
 * @param end_size The size in bytes of the end range.
 * @param is_empty The function sets it to `1` if the non-empty domain is
 *     empty (i.e., the array does not contain any data yet), and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_non_empty_domain_var_size_from_index(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint32_t idx,
    uint64_t* start_size,
    uint64_t* end_size,
    int32_t* is_empty) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain range sizes from an array for a given
 * dimension name. This is the union of the non-empty domains of the array
 * fragments on the given dimension. Applicable only to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t is_empty;
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 * uint64_t start_size, end_size;
 * tiledb_array_get_non_empty_domain_var_size_from_name(
 *     ctx, array, "d", &start_size, &end_size, &is_empty);
 * // If non-empty domain range is `[aa, dddd]`, then `start_size = 2`
 * // and `end_size = 4`.
 * @endcode
 *
 * @param ctx The TileDB context
 * @param array The array object (must be opened beforehand).
 * @param name The dimension name.
 * @param start_size The size in bytes of the start range.
 * @param end_size The size in bytes of the end range.
 * @param is_empty The function sets it to `1` if the non-empty domain is
 *     empty (i.e., the array does not contain any data yet), and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_non_empty_domain_var_size_from_name(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    const char* name,
    uint64_t* start_size,
    uint64_t* end_size,
    int32_t* is_empty) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain from an array for a given
 * dimension index. This is the union of the non-empty domains of the array
 * fragments on the given dimension. Applicable only to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t is_empty;
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 *
 * // Get range sizes first
 * uint64_t start_size, end_size;
 * tiledb_array_get_non_empty_domain_var_size_from_index(
 *     ctx, array, 0, &start_size, &end_size, &is_empty);
 *
 * // Get domain
 * char start[start_size];
 * char end[end_size];
 * tiledb_array_get_non_empty_domain_var_from_index(
 *     ctx, array, 0, start, end, &is_empty);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param array The array object (must be opened beforehand).
 * @param idx The dimension index, following the order as it was defined
 *      in the domain of the array schema.
 * @param start The domain range start to set.
 * @param end The domain range end to set.
 * @param is_empty The function sets it to `1` if the non-empty domain is
 *     empty (i.e., the array does not contain any data yet), and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_non_empty_domain_var_from_index(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint32_t idx,
    void* start,
    void* end,
    int32_t* is_empty) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain from an array for a given
 * dimension name. This is the union of the non-empty domains of the array
 * fragments on the given dimension. Applicable only to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t domain[2];
 * int32_t is_empty;
 * tiledb_array_t* array;
 * tiledb_array_alloc(ctx, "my_array", &array);
 * tiledb_array_open(ctx, array, TILEDB_READ);
 *
 * // Get range sizes first
 * uint64_t start_size, end_size;
 * tiledb_array_get_non_empty_domain_var_size_from_name(
 *     ctx, array, "d", &start_size, &end_size, &is_empty);
 *
 * // Get domain
 * char start[start_size];
 * char end[end_size];
 * tiledb_array_get_non_empty_domain_var_from_name(
 *     ctx, array, "d", start, end, &is_empty);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param array The array object (must be opened beforehand).
 * @param name The dimension name.
 * @param start The domain range start to set.
 * @param end The domain range end to set.
 * @param is_empty The function sets it to `1` if the non-empty domain is
 *     empty (i.e., the array does not contain any data yet), and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_non_empty_domain_var_from_name(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    const char* name,
    void* start,
    void* end,
    int32_t* is_empty) TILEDB_NOEXCEPT;

/**
 * Retrieves the URI the array was opened with. It outputs an error
 * if the array is not open.
 *
 * @param ctx The TileDB context.
 * @param array The input array.
 * @param array_uri The array URI to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_uri(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    const char** array_uri) TILEDB_NOEXCEPT;

/**
 * Retrieves the encryption type the array at the given URI was created with.
 *
 * @param ctx The TileDB context.
 * @param array_uri The array URI.
 * @param encryption_type The array encryption type to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_encryption_type(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_encryption_type_t* encryption_type) TILEDB_NOEXCEPT;

/**
 * It puts a metadata key-value item to an open array. The array must
 * be opened in WRITE mode, otherwise the function will error out.
 *
 * @param ctx The TileDB context.
 * @param array An array opened in WRITE mode.
 * @param key The key of the metadata item to be added. UTF-8 encodings
 *     are acceptable.
 * @param value_type The datatype of the value.
 * @param value_num The value may consist of more than one items of the
 *     same datatype. This argument indicates the number of items in the
 *     value component of the metadata.
 * @param value The metadata value in binary form.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note The writes will take effect only upon closing the array.
 */
TILEDB_EXPORT int32_t tiledb_array_put_metadata(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    const char* key,
    tiledb_datatype_t value_type,
    uint32_t value_num,
    const void* value) TILEDB_NOEXCEPT;

/**
 * It deletes a metadata key-value item from an open array. The array must
 * be opened in WRITE mode, otherwise the function will error out.
 *
 * @param ctx The TileDB context.
 * @param array An array opened in WRITE mode.
 * @param key The key of the metadata item to be deleted.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note The writes will take effect only upon closing the array.
 *
 * @note If the key does not exist, this will take no effect
 *     (i.e., the function will not error out).
 */
TILEDB_EXPORT int32_t tiledb_array_delete_metadata(
    tiledb_ctx_t* ctx, tiledb_array_t* array, const char* key) TILEDB_NOEXCEPT;

/**
 * It gets a metadata key-value item from an open array. The array must
 * be opened in READ mode, otherwise the function will error out.
 *
 * @param ctx The TileDB context.
 * @param array An array opened in READ mode.
 * @param key The key of the metadata item to be retrieved. UTF-8 encodings
 *     are acceptable.
 * @param value_type The datatype of the value.
 * @param value_num The value may consist of more than one items of the
 *     same datatype. This argument indicates the number of items in the
 *     value component of the metadata. Keys with empty values are indicated
 *     by value_num == 1 and value == NULL.
 * @param value The metadata value in binary form.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note If the key does not exist, then `value` will be NULL.
 */
TILEDB_EXPORT int32_t tiledb_array_get_metadata(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    const char* key,
    tiledb_datatype_t* value_type,
    uint32_t* value_num,
    const void** value) TILEDB_NOEXCEPT;

/**
 * It gets then number of metadata items in an open array. The array must
 * be opened in READ mode, otherwise the function will error out.
 *
 * @param ctx The TileDB context.
 * @param array An array opened in READ mode.
 * @param num The number of metadata items to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_metadata_num(
    tiledb_ctx_t* ctx, tiledb_array_t* array, uint64_t* num) TILEDB_NOEXCEPT;

/**
 * It gets a metadata item from an open array using an index.
 * The array must be opened in READ mode, otherwise the function will
 * error out.
 *
 * @param ctx The TileDB context.
 * @param array An array opened in READ mode.
 * @param index The index used to get the metadata.
 * @param key The metadata key.
 * @param key_len The metadata key length.
 * @param value_type The datatype of the value.
 * @param value_num The value may consist of more than one items of the
 *     same datatype. This argument indicates the number of items in the
 *     value component of the metadata.
 * @param value The metadata value in binary form.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_array_get_metadata_from_index(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    uint64_t index,
    const char** key,
    uint32_t* key_len,
    tiledb_datatype_t* value_type,
    uint32_t* value_num,
    const void** value) TILEDB_NOEXCEPT;

/**
 * Checks whether a key exists in metadata from an open array. The array must
 * be opened in READ mode, otherwise the function will error out.
 *
 * @param ctx The TileDB context.
 * @param array An array opened in READ mode.
 * @param key The key to be checked. UTF-8 encoding are acceptable.
 * @param value_type The datatype of the value, if any.
 * @param has_key Set to `1` if the metadata with given key exists, else `0`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note If the key does not exist, then `value` will be NULL.
 */
TILEDB_EXPORT int32_t tiledb_array_has_metadata_key(
    tiledb_ctx_t* ctx,
    tiledb_array_t* array,
    const char* key,
    tiledb_datatype_t* value_type,
    int32_t* has_key) TILEDB_NOEXCEPT;

/**
 * Consolidates the array metadata into a single array metadata file.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_consolidate_metadata(
 *     ctx, "hdfs:///tiledb_arrays/my_array", nullptr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The name of the TileDB array whose metadata will
 *     be consolidated.
 * @param config Configuration parameters for the consolidation
 *     (`nullptr` means default, which will use the config from `ctx`).
 * @return `TILEDB_OK` on success, and `TILEDB_ERR` on error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_consolidate_metadata(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_config_t* config) TILEDB_NOEXCEPT;

/**
 * Consolidates the array metadata of an encrypted array into a single file.
 *
 * **Example:**
 *
 * @code{.c}
 * uint8_t key[32] = ...;
 * tiledb_array_consolidate_metadata_with_key(
 *     ctx, "hdfs:///tiledb_arrays/my_array",
 *     TILEDB_AES_256_GCM, key, sizeof(key), nullptr);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The name of the TileDB array to be consolidated.
 * @param encryption_type The encryption type to use.
 * @param encryption_key The encryption key to use.
 * @param key_length Length in bytes of the encryption key.
 * @param config Configuration parameters for the consolidation
 *     (`nullptr` means default, which will use the config from `ctx`).
 *
 * @return `TILEDB_OK` on success, and `TILEDB_ERR` on error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_array_consolidate_metadata_with_key(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_encryption_type_t encryption_type,
    const void* encryption_key,
    uint32_t key_length,
    tiledb_config_t* config) TILEDB_NOEXCEPT;

/* ********************************* */
/*          OBJECT MANAGEMENT        */
/* ********************************* */

/**
 * Returns the TileDB object type for a given resource path.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_object_t type;
 * tiledb_object_type(ctx, "arrays/my_array", &type);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param path The URI path to the TileDB resource.
 * @param type The type to be retrieved.
 * @return `TILEDB_OK` on success, `TILEDB_ERR` on error.
 */
TILEDB_EXPORT int32_t tiledb_object_type(
    tiledb_ctx_t* ctx, const char* path, tiledb_object_t* type) TILEDB_NOEXCEPT;

/**
 * Deletes a TileDB resource (group, array, key-value).
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_object_remove(ctx, "arrays/my_array");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param path The URI path to the tiledb resource.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_object_remove(tiledb_ctx_t* ctx, const char* path)
    TILEDB_NOEXCEPT;

/**
 * Moves a TileDB resource (group, array, key-value).
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_object_move(ctx, "arrays/my_array", "arrays/my_array_2");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param old_path The old TileDB directory.
 * @param new_path The new TileDB directory.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_object_move(
    tiledb_ctx_t* ctx,
    const char* old_path,
    const char* new_path) TILEDB_NOEXCEPT;

/**
 * Walks (iterates) over the TileDB objects contained in *path*. The traversal
 * is done recursively in the order defined by the user. The user provides
 * a callback function which is applied on each of the visited TileDB objects.
 * The iteration continues for as long the callback returns non-zero, and stops
 * when the callback returns 0. Note that this function ignores any object
 * (e.g., file or directory) that is not TileDB-related.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_object_walk(ctx, "arrays", TILEDB_PREORDER, NULL, NULL);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param path The path in which the traversal will occur.
 * @param order The order of the recursive traversal (e.g., pre-order or
 *     post-order.
 * @param callback The callback function to be applied on every visited object.
 *     The callback should return `0` if the iteration must stop, and `1`
 *     if the iteration must continue. It takes as input the currently visited
 *     path, the type of that path (e.g., array or group), and the data
 *     provided by the user for the callback. The callback returns `-1` upon
 *     error. Note that `path` in the callback will be an **absolute** path.
 * @param data The data passed in the callback as the last argument.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_object_walk(
    tiledb_ctx_t* ctx,
    const char* path,
    tiledb_walk_order_t order,
    int32_t (*callback)(const char*, tiledb_object_t, void*),
    void* data) TILEDB_NOEXCEPT;

/**
 * Similar to `tiledb_walk`, but now the function visits only the children of
 * `path` (i.e., it does not recursively continue to the children directories).
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_object_ls(ctx, "arrays", NULL, NULL);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param path The path in which the traversal will occur.
 * @param callback The callback function to be applied on every visited object.
 *     The callback should return `0` if the iteration must stop, and `1`
 *     if the iteration must continue. It takes as input the currently visited
 *     path, the type of that path (e.g., array or group), and the data
 *     provided by the user for the callback. The callback returns `-1` upon
 *     error. Note that `path` in the callback will be an **absolute** path.
 * @param data The data passed in the callback as the last argument.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_object_ls(
    tiledb_ctx_t* ctx,
    const char* path,
    int32_t (*callback)(const char*, tiledb_object_t, void*),
    void* data) TILEDB_NOEXCEPT;

/* ****************************** */
/*        VIRTUAL FILESYSTEM      */
/* ****************************** */

/**
 * Creates a virtual filesystem object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_t* vfs;
 * tiledb_vfs_alloc(ctx, config, &vfs);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object to be created.
 * @param config Configuration parameters.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_alloc(
    tiledb_ctx_t* ctx,
    tiledb_config_t* config,
    tiledb_vfs_t** vfs) TILEDB_NOEXCEPT;

/**
 * Frees a virtual filesystem object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_free(&vfs);
 * @endcode
 *
 * @param vfs The virtual filesystem object to be freed.
 */
TILEDB_EXPORT void tiledb_vfs_free(tiledb_vfs_t** vfs) TILEDB_NOEXCEPT;

/**
 * Retrieves the config from a VFS context.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_config_t* config;
 * tiledb_vfs_get_config(ctx, vfs, &config);
 * // Make sure to free the retrieved config
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The VFS object.
 * @param config The config to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_get_config(
    tiledb_ctx_t* ctx,
    tiledb_vfs_t* vfs,
    tiledb_config_t** config) TILEDB_NOEXCEPT;

/**
 * Creates an object-store bucket.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_create_bucket(ctx, vfs, "s3://tiledb");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the bucket to be created.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_create_bucket(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri) TILEDB_NOEXCEPT;

/**
 * Deletes an object-store bucket.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_delete_bucket(ctx, vfs, "s3://tiledb");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the bucket to be deleted.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_remove_bucket(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri) TILEDB_NOEXCEPT;

/**
 * Deletes the contents of an object-store bucket.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_empty_bucket(ctx, vfs, "s3://tiledb");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the bucket to be emptied.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_empty_bucket(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri) TILEDB_NOEXCEPT;

/**
 * Checks if an object-store bucket is empty.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t is_empty;
 * tiledb_vfs_is_empty_bucket(ctx, vfs, "s3://tiledb", &empty);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the bucket.
 * @param is_empty Sets it to `1` if the input bucket is empty,
 *     and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_is_empty_bucket(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri, int32_t* is_empty)
    TILEDB_NOEXCEPT;

/**
 * Checks if an object-store bucket exists.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t exists;
 * tiledb_vfs_is_bucket(ctx, vfs, "s3://tiledb", &exists);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the bucket.
 * @param is_bucket Sets it to `1` if the input URI is a bucket, and `0`
 *     otherwise.
 * @return TILEDB_OK for success and TILEDB_ERR for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_is_bucket(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri, int32_t* is_bucket)
    TILEDB_NOEXCEPT;

/**
 * Creates a directory.
 *
 * - On S3, this is a noop.
 * - On all other backends, if the directory exists, the function
 *   just succeeds without doing anything.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_create_dir(ctx, vfs, "hdfs:///temp/my_dir");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the directory to be created.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_create_dir(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri) TILEDB_NOEXCEPT;

/**
 * Checks if a directory exists.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t exists;
 * tiledb_vfs_is_dir(ctx, vfs, "hdfs:///temp/my_dir", &exists);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the directory.
 * @param is_dir Sets it to `1` if the directory exists and `0`
 *     otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note For S3, this function will return `true` if there is an object
 *     with prefix `uri/` (TileDB will append `/` internally to `uri`
 *     only if it does not exist), and `false` othewise.
 */
TILEDB_EXPORT int32_t tiledb_vfs_is_dir(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri, int32_t* is_dir)
    TILEDB_NOEXCEPT;

/**
 * Removes a directory (recursively).
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_remove_dir(ctx, vfs, "hdfs:///temp/my_dir");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The uri of the directory to be removed
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_remove_dir(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri) TILEDB_NOEXCEPT;

/**
 * Checks if a file exists.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t exists;
 * tiledb_vfs_is_file(ctx, vfs, "hdfs:///temp/my_file", &is_file);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the file.
 * @param is_file Sets it to `1` if the file exists and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_is_file(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri, int32_t* is_file)
    TILEDB_NOEXCEPT;

/**
 * Deletes a file.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_remove_file(ctx, vfs, "hdfs:///temp/my_file");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the file.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_remove_file(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri) TILEDB_NOEXCEPT;

/**
 * Retrieves the size of a directory. This function is **recursive**, i.e.,
 * it will consider all files in the directory tree rooted at `uri`.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t dir_size;
 * tiledb_vfs_dir_size(ctx, vfs, "hdfs:///temp/my_dir", &dir_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the directory.
 * @param size The directory size to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_dir_size(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri, uint64_t* size)
    TILEDB_NOEXCEPT;

/**
 * Retrieves the size of a file.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t file_size;
 * tiledb_vfs_file_size(ctx, vfs, "hdfs:///temp/my_file", &file_size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the file.
 * @param size The file size to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_file_size(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri, uint64_t* size)
    TILEDB_NOEXCEPT;

/**
 * Renames a file. If the destination file exists, it will be overwritten.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_move_file(
 * ctx, vfs, "hdfs:///temp/my_file", "hdfs::///new_file");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param old_uri The old URI.
 * @param new_uri The new URI.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_move_file(
    tiledb_ctx_t* ctx,
    tiledb_vfs_t* vfs,
    const char* old_uri,
    const char* new_uri) TILEDB_NOEXCEPT;

/**
 * Renames a directory.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_move_dir(ctx, vfs, "hdfs:///temp/my_dir", "hdfs::///new_dir");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param old_uri The old URI.
 * @param new_uri The new URI.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_move_dir(
    tiledb_ctx_t* ctx,
    tiledb_vfs_t* vfs,
    const char* old_uri,
    const char* new_uri) TILEDB_NOEXCEPT;

/**
 * Copies a file. If the destination file exists, it will be overwritten.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_copy_file(
 * ctx, vfs, "hdfs:///temp/my_file", "hdfs::///new_file");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param old_uri The old URI.
 * @param new_uri The new URI.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_copy_file(
    tiledb_ctx_t* ctx,
    tiledb_vfs_t* vfs,
    const char* old_uri,
    const char* new_uri) TILEDB_NOEXCEPT;

/**
 * Copies a directory. If the destination directory exists, it will be
 * overwritten.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_copy_dir(
 *  ctx, vfs, "hdfs:///temp/my_dir", "hdfs::///new_dir");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param old_uri The old URI.
 * @param new_uri The new URI.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_copy_dir(
    tiledb_ctx_t* ctx,
    tiledb_vfs_t* vfs,
    const char* old_uri,
    const char* new_uri) TILEDB_NOEXCEPT;

/**
 * Prepares a file for reading/writing.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_fh_t* fh;
 * tiledb_vfs_open(ctx, vfs, "some_file", TILEDB_VFS_READ, &fh);
 * // Make sure to close and delete the created file handle
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the file.
 * @param mode The mode in which the file is opened:
 *     - `TILEDB_VFS_READ`:
 *       The file is opened for reading. An error is returned if the file
 *       does not exist.
 *     - `TILEDB_VFS_WRITE`:
 *       The file is opened for writing. If the file exists, it will be
 *       overwritten.
 *     - `TILEDB_VFS_APPEND`:
 *       The file is opened for writing. If the file exists, the write
 *       will start from the end of the file. Note that S3 does not
 *       support this operation and, thus, an error will be thrown in
 *       that case.
 * @param fh The file handle that is created. This will be used in
 *     `tiledb_vfs_read`, `tiledb_vfs_write` and `tiledb_vfs_sync`.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` or `TILEDB_OOM` for error.
 *
 * @note If the file is closed after being opened, without having
 *     written any data to it, the file will not be created. If you
 *     wish to create an empty file, use `tiledb_vfs_touch`
 *     instead.
 */
TILEDB_EXPORT int32_t tiledb_vfs_open(
    tiledb_ctx_t* ctx,
    tiledb_vfs_t* vfs,
    const char* uri,
    tiledb_vfs_mode_t mode,
    tiledb_vfs_fh_t** fh) TILEDB_NOEXCEPT;

/**
 * Closes a file. This is flushes the buffered data into the file
 * when the file was opened in write (or append) mode. It is particularly
 * important to be called after S3 writes, as otherwise the writes will
 * not take effect.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_close(ctx, vfs, fh);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fh The file handle.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_close(tiledb_ctx_t* ctx, tiledb_vfs_fh_t* fh)
    TILEDB_NOEXCEPT;

/**
 * Reads from a file.
 *
 * **Example:**
 *
 * @code{.c}
 * char buffer[10000];
 * tiledb_vfs_read(ctx, fh, 100, buffer, 10000);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fh The URI file handle.
 * @param offset The offset in the file where the read begins.
 * @param buffer The buffer to read into.
 * @param nbytes Number of bytes to read.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_read(
    tiledb_ctx_t* ctx,
    tiledb_vfs_fh_t* fh,
    uint64_t offset,
    void* buffer,
    uint64_t nbytes) TILEDB_NOEXCEPT;

/**
 * Writes the contents of a buffer into a file. Note that this
 * function only **appends** data at the end of the file. If the
 * file does not exist, it will be created.
 *
 * **Example:**
 *
 * @code{.c}
 * const char* msg = "This will be written to the file";
 * tiledb_vfs_write(ctx, fh, buffer, strlen(msg));
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fh The URI file handle.
 * @param buffer The buffer to write from.
 * @param nbytes Number of bytes to write.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_write(
    tiledb_ctx_t* ctx, tiledb_vfs_fh_t* fh, const void* buffer, uint64_t nbytes)
    TILEDB_NOEXCEPT;

/**
 * Syncs (flushes) a file.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_sync(ctx, fh);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fh The URI file handle.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note This has no effect for S3.
 */
TILEDB_EXPORT int32_t tiledb_vfs_sync(tiledb_ctx_t* ctx, tiledb_vfs_fh_t* fh)
    TILEDB_NOEXCEPT;

/**
 * The function visits only the children of `path` (i.e., it does not
 * recursively continue to the children directories) and applies the `callback`
 * function using the input `data`.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_ls(ctx, vfs, "my_dir", NULL, NULL);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param path The path in which the traversal will occur.
 * @param callback The callback function to be applied on every visited object.
 *     The callback should return `0` if the iteration must stop, and `1`
 *     if the iteration must continue. It takes as input the currently visited
 *     path, the type of that path (e.g., array or group), and the data
 *     provided by the user for the callback. The callback returns `-1` upon
 *     error. Note that `path` in the callback will be an **absolute** path.
 * @param data The data passed in the callback as the last argument.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_ls(
    tiledb_ctx_t* ctx,
    tiledb_vfs_t* vfs,
    const char* path,
    int32_t (*callback)(const char*, void*),
    void* data) TILEDB_NOEXCEPT;

/**
 * Frees a file handle.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_fh_free(&fh);
 * @endcode
 *
 * @param fh The URI file handle.
 */
TILEDB_EXPORT void tiledb_vfs_fh_free(tiledb_vfs_fh_t** fh) TILEDB_NOEXCEPT;

/**
 * Checks if a file handle is closed.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t is_closed;
 * tiledb_vfs_fh_is_closed(ctx, fh, &is_closed);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fh The URI file handle.
 * @param is_closed Set to `1` if the file handle is closed, and `0` otherwise.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_fh_is_closed(
    tiledb_ctx_t* ctx, tiledb_vfs_fh_t* fh, int32_t* is_closed) TILEDB_NOEXCEPT;

/**
 * Touches a file, i.e., creates a new empty file.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_vfs_touch(ctx, vfs, "my_empty_file");
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param vfs The virtual filesystem object.
 * @param uri The URI of the file to be created.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_vfs_touch(
    tiledb_ctx_t* ctx, tiledb_vfs_t* vfs, const char* uri) TILEDB_NOEXCEPT;

/* ****************************** */
/*              URI               */
/* ****************************** */

/**
 * Converts the given file URI to a null-terminated platform-native file path.
 *
 * **Example:**
 *
 * @code{.c}
 * char path[TILEDB_MAX_PATH];
 * uint32_t length = TILEDB_MAX_PATH; // Must be set to non-zero
 * tiledb_uri_to_path(ctx, "file:///my_array", path, &length);
 * // This will set "my_array" to `path`
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param uri The URI to be converted.
 * @param path_out The buffer where the converted path string is stored.
 * @param path_length The length of the path buffer. On return, this is set to
 * the length of the converted path string, or 0 on error.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 *
 * @note The path_out buffer must be allocated according to the platform's
 * maximum path length (e.g. `TILEDB_MAX_PATH), which includes space for the
 * terminating null character.
 */
TILEDB_EXPORT int32_t tiledb_uri_to_path(
    tiledb_ctx_t* ctx, const char* uri, char* path_out, uint32_t* path_length)
    TILEDB_NOEXCEPT;

/* ****************************** */
/*             Stats              */
/* ****************************** */

/**
 * Enable internal statistics gathering.
 *
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_stats_enable(void) TILEDB_NOEXCEPT;

/**
 * Disable internal statistics gathering.
 *
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_stats_disable(void) TILEDB_NOEXCEPT;

/**
 * Reset all internal statistics counters to 0.
 *
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_stats_reset(void) TILEDB_NOEXCEPT;

/**
 * Dump all internal statistics counters to some output (e.g.,
 * file or stdout).
 *
 * @param out The output.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_stats_dump(FILE* out) TILEDB_NOEXCEPT;

/**
 * Dump all internal statistics counters to an output string. The caller is
 * responsible for freeing the resulting string.
 *
 * **Example:**
 *
 * @code{.c}
 * char *stats_str;
 * tiledb_stats_dump_str(&stats_str);
 * // ...
 * tiledb_stats_free_str(&stats_str);
 * @endcode
 *
 * @param out Will be set to point to an allocated string containing the stats.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_stats_dump_str(char** out) TILEDB_NOEXCEPT;

/**
 * Dump all raw internal statistics counters to some output (e.g.,
 * file or stdout) as a JSON.
 *
 * @param out The output.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_stats_raw_dump(FILE* out) TILEDB_NOEXCEPT;

/**
 * Dump all raw internal statistics counters to a JSON-formatted output string.
 * The caller is responsible for freeing the resulting string.
 *
 * **Example:**
 *
 * @code{.c}
 * char *stats_str;
 * tiledb_stats_raw_dump_str(&stats_str);
 * // ...
 * tiledb_stats_raw_free_str(&stats_str);
 * @endcode
 *
 * @param out Will be set to point to an allocated string containing the stats.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_stats_raw_dump_str(char** out) TILEDB_NOEXCEPT;

/**
 *
 * Free the memory associated with a previously dumped stats string.
 *
 * @param out Pointer to a previously allocated stats string.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_stats_free_str(char** out) TILEDB_NOEXCEPT;

/* ****************************** */
/*          Heap Profiler         */
/* ****************************** */

/**
 * Enable heap profiling.
 *
 * @param file_name_prefix If empty or null, stats are dumped
 *   to stdout. If non-empty, this specifies the file_name prefix to
 *   write to. For example, value "tiledb_mem_stats" will write
 *   to "tiledb_mem_stats__1611170501", where the postfix is
 *   determined by the current epoch.
 * @param dump_interval_ms If non-zero, this spawns a dedicated
 *   thread to dump on this time interval.
 * @param dump_interval_bytes If non-zero, a dump will occur when
 *   the total number of lifetime allocated bytes is increased by
 *   more than this amount.
 * @param dump_threshold_bytes If non-zero, labeled allocations with
 *   a number of bytes lower than this threshold will not be reported
 *   in the dump.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_heap_profiler_enable(
    const char* file_name_prefix,
    uint64_t dump_interval_ms,
    uint64_t dump_interval_bytes,
    uint64_t dump_threshold_bytes) TILEDB_NOEXCEPT;

/* ****************************** */
/*          FRAGMENT INFO         */
/* ****************************** */

/**
 * Creates a fragment info object for a given array, and fetches all
 * the fragment information for that array.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_fragment_info* fragment_info;
 * tiledb_fragment_info_alloc(ctx, "array_uri", &fragment_info);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param array_uri The array URI.
 * @param fragment_info The fragment info object to be created and populated.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_alloc(
    tiledb_ctx_t* ctx,
    const char* array_uri,
    tiledb_fragment_info_t** fragment_info) TILEDB_NOEXCEPT;

/**
 * Frees a fragment info object.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_fragment_info_free(&fragment_info);
 * @endcode
 *
 * @param fragment_info The fragment info object to be freed.
 */
TILEDB_EXPORT void tiledb_fragment_info_free(
    tiledb_fragment_info_t** fragment_info) TILEDB_NOEXCEPT;

/**
 * Set the fragment info config. Useful for passing timestamp ranges and
 * encryption key via the config before loading the fragment info.
 *
 *  * **Example:**
 *
 * @code{.c}
 * tiledb_fragment_info* fragment_info;
 * tiledb_fragment_info_alloc(ctx, "array_uri", &fragment_info);
 *
 * tiledb_config_t* config;
 * tiledb_error_t* error = NULL;
 * tiledb_config_alloc(&config, &error);
 * tiledb_config_set(config, "sm.tile_cache_size", "1000000", &error);
 *
 * tiledb_fragment_info_load(ctx, fragment_info);
 * @endcode

 */
TILEDB_EXPORT int32_t tiledb_fragment_info_set_config(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    tiledb_config_t* config) TILEDB_NOEXCEPT;

/**
 * Retrieves the config from fragment info.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_config_t* config;
 * tiledb_fragment_info_get_config(ctx, vfs, &config);
 * // Make sure to free the retrieved config
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param config The config to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_OOM` or `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_config(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    tiledb_config_t** config) TILEDB_NOEXCEPT;

/**
 * Loads the fragment info.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_fragment_info_load(ctx, fragment_info);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_load(
    tiledb_ctx_t* ctx, tiledb_fragment_info_t* fragment_info) TILEDB_NOEXCEPT;

/**
 * Loads the fragment info from an encrypted array.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_fragment_info_load_with_key(
 *     ctx, fragment_info, TILEDB_AES_256_GCM, key, sizeof(key));
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param encryption_type The encryption type to use.
 * @param encryption_key The encryption key to use.
 * @param key_length Length in bytes of the encryption key.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_DEPRECATED_EXPORT int32_t tiledb_fragment_info_load_with_key(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    tiledb_encryption_type_t encryption_type,
    const void* encryption_key,
    uint32_t key_length) TILEDB_NOEXCEPT;

/**
 * Gets a fragment name.
 *
 * **Example:**
 *
 * @code{.c}
 * const char* name;
 * tiledb_fragment_info_get_fragment_name(ctx, fragment_info, 1, &name);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param name The fragment name to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_fragment_name(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    const char** name) TILEDB_NOEXCEPT;

/**
 * Gets the number of fragments.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t fragment_num;
 * tiledb_fragment_info_get_fragment_num(ctx, fragment_info, &fragment_num);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param fragment_num The number of fragments to retrieve.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_fragment_num(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t* fragment_num) TILEDB_NOEXCEPT;

/**
 * Gets a fragment URI.
 *
 * **Example:**
 *
 * @code{.c}
 * const char* uri;
 * tiledb_fragment_info_get_fragment_uri(ctx, fragment_info, 1, &uri);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param uri The fragment URI to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_fragment_uri(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    const char** uri) TILEDB_NOEXCEPT;

/**
 * Gets the fragment size in bytes.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t size;
 * tiledb_fragment_info_get_fragment_size(ctx, fragment_info, 1, &size);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param size The fragment size to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_fragment_size(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint64_t* size) TILEDB_NOEXCEPT;

/**
 * Checks if a fragment is dense.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t dense;
 * tiledb_fragment_info_get_dense(ctx, fragment_info, 1, &dense);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param dense `1` if the fragment is dense.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_dense(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    int32_t* dense) TILEDB_NOEXCEPT;

/**
 * Checks if a fragment is sparse.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t sparse;
 * tiledb_fragment_info_get_sparse(ctx, fragment_info, 1, &sparse);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param sparse `1` if the fragment is sparse.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_sparse(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    int32_t* sparse) TILEDB_NOEXCEPT;

/**
 * Gets the timestamp range of a fragment.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t start, end;
 * tiledb_fragment_info_get_timestamp_range(ctx, fragment_info, 1, &start,
 * &end);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param start The start of the timestamp range to be retrieved.
 * @param end The end of the timestamp range to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_timestamp_range(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint64_t* start,
    uint64_t* end) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain from a given fragment for a given
 * dimension index.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t domain[2];
 * tiledb_fragment_info_get_non_empty_domain_from_index(
 *     ctx, fragment_info, 0, 0, domain);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param did The dimension index, following the order as it was defined
 *      in the domain of the array schema.
 * @param domain The domain to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_non_empty_domain_from_index(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t did,
    void* domain) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain from a given fragment for a given
 * dimension name.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t domain[2];
 * tiledb_fragment_info_get_non_empty_domain_from_name(
 *     ctx, fragment_info, 0, "d1", domain);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param dim_name The dimension name.
 * @param domain The domain to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_non_empty_domain_from_name(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    const char* dim_name,
    void* domain) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain range sizes from a fragment for a given
 * dimension index. Applicable to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t start_size, end_size;
 * tiledb_fragment_info_get_non_empty_domain_var_size_from_index(
 *     ctx, fragment_info, 0, &start_size, &end_size);
 * // If non-empty domain range is `[aa, dddd]`, then `start_size = 2`
 * // and `end_size = 4`.
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment information object.
 * @param fid The fragment index.
 * @param did The dimension index, following the order as it was defined
 *      in the domain of the array schema.
 * @param start_size The size in bytes of the start range.
 * @param end_size The size in bytes of the end range.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t
tiledb_fragment_info_get_non_empty_domain_var_size_from_index(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t did,
    uint64_t* start_size,
    uint64_t* end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain range sizes from a fragment for a given
 * dimension name. Applicable to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t start_size, end_size;
 * tiledb_fragment_info_get_non_empty_domain_var_size_from_name(
 *     ctx, fragment_info, "d", &start_size, &end_size);
 * // If non-empty domain range is `[aa, dddd]`, then `start_size = 2`
 * // and `end_size = 4`.
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment information object.
 * @param fid The fragment index.
 * @param dim_name The dimension name.
 * @param start_size The size in bytes of the start range.
 * @param end_size The size in bytes of the end range.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t
tiledb_fragment_info_get_non_empty_domain_var_size_from_name(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    const char* dim_name,
    uint64_t* start_size,
    uint64_t* end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain from a fragment for a given
 * dimension index. Applicable to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 *
 * // Get range sizes first
 * uint64_t start_size, end_size;
 * tiledb_fragment_info_get_non_empty_domain_var_size_from_index(
 *     ctx, fragment_info, 0, 0, &start_size, &end_size);
 *
 * // Get domain
 * char start[start_size];
 * char end[end_size];
 * tiledb_fragment_info_get_non_empty_domain_var_from_index(
 *     ctx, fragment_info, 0, 0, start, end);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The fragment index.
 * @param did The dimension index, following the order as it was defined
 *      in the domain of the array schema.
 * @param start The domain range start to set.
 * @param end The domain range end to set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_non_empty_domain_var_from_index(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t did,
    void* start,
    void* end) TILEDB_NOEXCEPT;

/**
 * Retrieves the non-empty domain from a fragment for a given dimension name.
 * Applicable to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 *
 * // Get range sizes first
 * uint64_t start_size, end_size;
 * tiledb_fragment_info_get_non_empty_domain_var_size_from_name(
 *     ctx, fragment_info, 0, "d", &start_size, &end_size);
 *
 * // Get domain
 * char start[start_size];
 * char end[end_size];
 * tiledb_fragment_info_get_non_empty_domain_var_from_name(
 *     ctx, fragment_info, 0, "d", start, end);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The fragment index.
 * @param dim_name The dimension name.
 * @param start The domain range start to set.
 * @param end The domain range end to set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_non_empty_domain_var_from_name(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    const char* dim_name,
    void* start,
    void* end) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of MBRs from the fragment.
 *
 * In the case of sparse fragments, this is the number of physical tiles.
 *
 * Dense fragments do not contain MBRs.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t mbr_num;
 * tiledb_fragment_info_get_mbr_num(ctx, fragment_info, 0, &mbr_num);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param mbrs_num The number of MBRs to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_mbr_num(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint64_t* mbr_num) TILEDB_NOEXCEPT;

/**
 * Retrieves the MBR from a given fragment for a given dimension index.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t mbr[2];
 * tiledb_fragment_info_get_mbr_from_index(ctx, fragment_info, 0, 0, 0, mbr);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param mid The mbr of the fragment of interest.
 * @param did The dimension index, following the order as it was defined
 *      in the domain of the array schema.
 * @param mbr The mbr to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_mbr_from_index(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t mid,
    uint32_t did,
    void* mbr) TILEDB_NOEXCEPT;

/**
 * Retrieves the MBR from a given fragment for a given dimension name.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t mbr[2];
 * tiledb_fragment_info_get_mbr_from_name(ctx, fragment_info, 0, 0, "d1", mbr);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param mid The mbr of the fragment of interest.
 * @param dim_name The dimension name.
 * @param mbr The mbr to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_mbr_from_name(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t mid,
    const char* dim_name,
    void* mbr) TILEDB_NOEXCEPT;

/**
 * Retrieves the MBR sizes from a fragment for a given dimension index.
 * Applicable to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t start_size, end_size;
 * tiledb_fragment_info_get_mbr_var_size_from_index(
 *     ctx, fragment_info, 0, 0, 0, &start_size, &end_size);
 * // If non-empty domain range is `[aa, dddd]`, then `start_size = 2`
 * // and `end_size = 4`.
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment information object.
 * @param fid The fragment index.
   @param mid The mbr of the fragment of interest.
 * @param did The dimension index, following the order as it was defined
 *      in the domain of the array schema.
 * @param start_size The size in bytes of the start range.
 * @param end_size The size in bytes of the end range.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_mbr_var_size_from_index(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t mid,
    uint32_t did,
    uint64_t* start_size,
    uint64_t* end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves the MBR range sizes from a fragment for a given dimension name.
 * Applicable to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t start_size, end_size;
 * tiledb_fragment_info_get_mbr_var_size_from_name(
 *     ctx, fragment_info, 0, 0, "d1", &start_size, &end_size);
 * // If non-empty domain range is `[aa, dddd]`, then `start_size = 2`
 * // and `end_size = 4`.
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment information object.
 * @param fid The fragment index.
 * @param mid The mbr of the fragment of interest.
 * @param dim_name The dimension name.
 * @param start_size The size in bytes of the start range.
 * @param end_size The size in bytes of the end range.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_mbr_var_size_from_name(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t mid,
    const char* dim_name,
    uint64_t* start_size,
    uint64_t* end_size) TILEDB_NOEXCEPT;

/**
 * Retrieves the MBR from a fragment for a given dimension index.
 * Applicable to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 *
 * // Get range sizes first
 * uint64_t start_size, end_size;
 * tiledb_fragment_info_get_mbr_var_size_from_index(
 *     ctx, fragment_info, 0, 0, 0, &start_size, &end_size);
 *
 * // Get domain
 * char start[start_size];
 * char end[end_size];
 * tiledb_fragment_info_get_mbr_var_from_index(
 *     ctx, fragment_info, 0, 0, 0, start, end);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The fragment index.
 * @param mid The mbr of the fragment of interest.
 * @param did The dimension index, following the order as it was defined
 *      in the domain of the array schema.
 * @param start The domain range start to set.
 * @param end The domain range end to set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_mbr_var_from_index(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t mid,
    uint32_t did,
    void* start,
    void* end) TILEDB_NOEXCEPT;

/**
 * Retrieves the MBR from a fragment for a given dimension name.
 * Applicable to var-sized dimensions.
 *
 * **Example:**
 *
 * @code{.c}
 *
 * // Get range sizes first
 * uint64_t start_size, end_size;
 * tiledb_fragment_info_get_mbr_var_size_from_name(
 *     ctx, fragment_info, 0, 0, "d1", &start_size, &end_size);
 *
 * // Get domain
 * char start[start_size];
 * char end[end_size];
 * tiledb_fragment_info_get_mbr_var_from_name(
 *     ctx, fragment_info, 0, 0, "d1", start, end);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The fragment index.
 * @param mid The mbr of the fragment of interest.
 * @param dim_name The dimension name.
 * @param start The domain range start to set.
 * @param end The domain range end to set.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_mbr_var_from_name(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t mid,
    const char* dim_name,
    void* start,
    void* end) TILEDB_NOEXCEPT;

/**
 * Retrieves the number of cells written to the fragment by the user.
 *
 * In the case of sparse fragments, this is the number of non-empty
 * cells in the fragment.
 *
 * In the case of dense fragments, TileDB may add fill
 * values to populate partially populated tiles. Those fill values
 * are counted in the returned number of cells. In other words,
 * the cell number is derived from the number of *integral* tiles
 * written in the file.
 *
 * **Example:**
 *
 * @code{.c}
 * uint64_t cell_num;
 * tiledb_fragment_info_get_cell_num(ctx, fragment_info, 0, &cell_num);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param cell_num The number of cells to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_cell_num(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint64_t* cell_num) TILEDB_NOEXCEPT;

/**
 * Retrieves the format version of a fragment.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t version;
 * tiledb_fragment_info_get_version(ctx, fragment_info, 0, &version);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param version The format version to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_version(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    uint32_t* version) TILEDB_NOEXCEPT;

/**
 * Checks if a fragment has consolidated metadata.
 *
 * **Example:**
 *
 * @code{.c}
 * int32_t has;
 * tiledb_fragment_info_has_consolidated_metadata(ctx, fragment_info, 0, &has);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param has True if the fragment has consolidated metadata.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_has_consolidated_metadata(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    int32_t* has) TILEDB_NOEXCEPT;

/**
 * Gets the number of fragments with unconsolidated metadata.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t unconsolidated;
 * tiledb_fragment_info_get_unconsolidated_metadata_num(ctx, fragment_info,
 * &unconsolidated);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param unconsolidated The number of fragments with unconsolidated metadata.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_unconsolidated_metadata_num(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t* unconsolidated) TILEDB_NOEXCEPT;

/**
 * Gets the number of fragments to vacuum.
 *
 * **Example:**
 *
 * @code{.c}
 * uint32_t to_vacuum_num;
 * tiledb_fragment_info_get_to_vacuum_num(ctx, fragment_info, &to_vacuum_num);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param to_vacuum_num The number of fragments to vacuum.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_to_vacuum_num(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t* to_vacuum_num) TILEDB_NOEXCEPT;

/**
 * Gets the URI of the fragment to vacuum with the given index.
 *
 * **Example:**
 *
 * @code{.c}
 * const char* uri;
 * tiledb_fragment_info_get_to_vacuum_uri(ctx, fragment_info, 1, &uri);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment to vacuum of interest.
 * @param uri The fragment URI to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_to_vacuum_uri(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    const char** uri) TILEDB_NOEXCEPT;

/**
 * Retrieves the array schema name a fragment.
 *
 * **Example:**
 *
 * @code{.c}
 * tiledb_array_schema_t* array_schema;
 * tiledb_fragment_info_get_array_schema(ctx, fragment_info, 0, &array_schema);
 * @endcode
 *
 * @param ctx The TileDB context
 * @param fragment_info The fragment info object.
 * @param fid The index of the fragment of interest.
 * @param array_schema The array schema to be retrieved.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_array_schema(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    tiledb_array_schema_t** array_schema) TILEDB_NOEXCEPT;

/**
 * Get the fragment info schema name.
 *
 * **Example:**
 *
 * @code{.c}
 * char* name;
 * tiledb_fragment_info_schema_name(ctx, fragment_info, &schema_name);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param name The schema name.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_get_array_schema_name(
    tiledb_ctx_t* ctx,
    tiledb_fragment_info_t* fragment_info,
    uint32_t fid,
    const char** schema_name) TILEDB_NOEXCEPT;

/**
 * Dumps the fragment info in ASCII format in the selected output.
 *
 * **Example:**
 *
 * The following prints the fragment info dump in standard output.
 *
 * @code{.c}
 * tiledb_fragment_info_dump(ctx, fragment_info, stdout);
 * @endcode
 *
 * @param ctx The TileDB context.
 * @param fragment_info The fragment info object.
 * @param out The output.
 * @return `TILEDB_OK` for success and `TILEDB_ERR` for error.
 */
TILEDB_EXPORT int32_t tiledb_fragment_info_dump(
    tiledb_ctx_t* ctx,
    const tiledb_fragment_info_t* fragment_info,
    FILE* out) TILEDB_NOEXCEPT;

#ifdef __cplusplus
}
#endif

#endif  // TILEDB_H
