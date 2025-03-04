/**
 * @file   vfs.h
 *
 * @section LICENSE
 *
 * The MIT License
 *
 * @copyright Copyright (c) 2017-2022 TileDB, Inc.
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
 * This file declares the VFS class.
 */

#ifndef TILEDB_VFS_H
#define TILEDB_VFS_H

#include <functional>
#include <list>
#include <set>
#include <string>
#include <vector>

#include "tiledb/common/common.h"
#include "tiledb/common/filesystem/directory_entry.h"
#include "tiledb/common/macros.h"
#include "tiledb/common/status.h"
#include "tiledb/sm/buffer/buffer.h"
#include "tiledb/sm/cache/lru_cache.h"
#include "tiledb/sm/config/config.h"
#include "tiledb/sm/filesystem/mem_filesystem.h"
#include "tiledb/sm/misc/cancelable_tasks.h"
#include "tiledb/sm/stats/stats.h"
#include "uri.h"

#ifdef _WIN32
#include "tiledb/sm/filesystem/win.h"
#else
#include "tiledb/sm/filesystem/posix.h"
#endif

#ifdef HAVE_GCS
#include "tiledb/sm/filesystem/gcs.h"
#endif

#ifdef HAVE_S3
#include "tiledb/sm/filesystem/s3.h"
#endif

#ifdef HAVE_HDFS
#include "tiledb/sm/filesystem/hdfs_filesystem.h"
#endif

#ifdef HAVE_AZURE
#include "tiledb/sm/filesystem/azure.h"
#endif

using namespace tiledb::common;

namespace tiledb::sm {

class Tile;

enum class Filesystem : uint8_t;
enum class VFSMode : uint8_t;

/** The VFS configuration parameters. */
struct VFSParameters {
  VFSParameters() = delete;

  VFSParameters(const Config& config)
      : max_batch_size_(config.get<uint64_t>("vfs.max_batch_size").value())
      , min_batch_gap_(config.get<uint64_t>("vfs.min_batch_gap").value())
      , min_batch_size_(config.get<uint64_t>("vfs.min_batch_size").value())
      , min_parallel_size_(
            config.get<uint64_t>("vfs.min_parallel_size").value())
      , read_ahead_cache_size_(
            config.get<uint64_t>("vfs.read_ahead_cache_size").value())
      , read_ahead_size_(config.get<uint64_t>("vfs.read_ahead_size").value()){};

  ~VFSParameters() = default;

  /** The maximum number of bytes in a batched read operation. */
  uint64_t max_batch_size_;

  /** The minimum number of bytes between two read batches. */
  uint64_t min_batch_gap_;

  /** The minimum number of bytes in a batched read operation. */
  uint64_t min_batch_size_;

  /** The minimum number of bytes in a parallel operation. */
  uint64_t min_parallel_size_;

  /** The byte size of the read-ahead cache. */
  uint64_t read_ahead_cache_size_;

  /** The byte size to read-ahead for each read. */
  uint64_t read_ahead_size_;
};

/**
 * This class implements a virtual filesystem that directs filesystem-related
 * function execution to the appropriate backend based on the input URI.
 */
class VFS {
 public:
  /* ********************************* */
  /*          TYPE DEFINITIONS         */
  /* ********************************* */
  /**
   * Multipart upload state definition used in the serialization of remote
   * global order writes. This state is a generalization of
   * the multipart upload state types currently defined independently by each
   * backend implementation.
   */
  struct MultiPartUploadState {
    struct CompletedParts {
      optional<std::string> e_tag;
      uint64_t part_number;
    };

    uint64_t part_number;
    optional<std::string> upload_id;
    std::vector<CompletedParts> completed_parts;
    Status status;
  };

  /* ********************************* */
  /*     CONSTRUCTORS & DESTRUCTORS    */
  /* ********************************* */

  /** Constructor.
   * @param parent_stats The parent stats to inherit from.
   * @param compute_tp Thread pool for compute-bound tasks.
   * @param io_tp Thread pool for io-bound tasks.
   * @param config Configuration parameters.
   **/
  VFS(stats::Stats* parent_stats,
      ThreadPool* compute_tp,
      ThreadPool* io_tp,
      const Config& config);

  /** Destructor. */
  ~VFS() = default;

  DISABLE_COPY_AND_COPY_ASSIGN(VFS);
  DISABLE_MOVE_AND_MOVE_ASSIGN(VFS);

  /* ********************************* */
  /*               API                 */
  /* ********************************* */

  /**
   * Returns the absolute path of the input string (mainly useful for
   * posix URI's).
   *
   * @param path The input path.
   * @return The string with the absolute path.
   */
  static std::string abs_path(const std::string& path);

  /**
   * Return a config object containing the VFS parameters. All other non-VFS
   * parameters will are set to default values.
   */
  Config config() const;

  /**
   * Creates a directory.
   *
   * - On S3, this is a noop.
   * - On all other backends, if the directory exists, the function
   *   just succeeds without doing anything.
   *
   * @param uri The URI of the directory.
   * @return Status
   */
  Status create_dir(const URI& uri) const;

  /**
   * Creates an empty file.
   *
   * @param uri The URI of the file.
   * @return Status
   */
  Status touch(const URI& uri) const;

  /**
   * Cancels all background or queued tasks.
   */
  Status cancel_all_tasks();

  /**
   * Creates an object store bucket.
   *
   * @param uri The name of the bucket to be created.
   * @return Status
   */
  Status create_bucket(const URI& uri) const;

  /**
   * Returns the size of the files in the input directory.
   * This function is **recursive**, i.e., it will calculate
   * the sum of the files in the entire directory tree rooted
   * at `dir_name`.
   *
   * @param dir_name The input directory.
   * @param dir_size The directory size to be retrieved, as the
   *     sum of the files one level deep.
   * @return Status
   */
  Status dir_size(const URI& dir_name, uint64_t* dir_size) const;

  /**
   * Deletes an object store bucket.
   *
   * @param uri The name of the bucket to be deleted.
   * @return Status
   */
  Status remove_bucket(const URI& uri) const;

  /**
   * Deletes the contents of an object store bucket.
   *
   * @param uri The name of the bucket to be emptied.
   * @return Status
   */
  Status empty_bucket(const URI& uri) const;

  /**
   * Removes a given directory (recursive)
   *
   * @param uri The uri of the directory to be removed
   * @return Status
   */
  Status remove_dir(const URI& uri) const;

  /**
   * Deletes a file.
   *
   * @param uri The URI of the file.
   * @return Status
   */
  Status remove_file(const URI& uri) const;

  /**
   * Retrieves the size of a file.
   *
   * @param uri The URI of the file.
   * @param size The file size to be retrieved.
   * @return Status
   */
  Status file_size(const URI& uri, uint64_t* size) const;

  /**
   * Checks if a directory exists.
   *
   * @param uri The URI of the directory.
   * @param is_dir Set to `true` if the directory exists and `false` otherwise.
   * @return Status
   *
   * @note For S3, this function will return `true` if there is an object
   *     with prefix `uri/` (TileDB will append `/` internally to `uri`
   *     only if it does not exist), and `false` othewise.
   */
  Status is_dir(const URI& uri, bool* is_dir) const;

  /**
   * Checks if a file exists.
   *
   * @param uri The URI of the file.
   * @param is_file Set to `true` if the file exists and `false` otherwise.
   * @return Status
   */
  Status is_file(const URI& uri, bool* is_file) const;

  /**
   * Checks if an object store bucket exists.
   *
   * @param uri The name of the object store bucket.
   * @return is_bucket Set to `true` if the bucket exists and `false` otherwise.
   * @return Status
   */
  Status is_bucket(const URI& uri, bool* is_bucket) const;

  /**
   * Checks if an object-store bucket is empty.
   *
   * @param uri The name of the object store bucket.
   * @param is_empty Set to `true` if the bucket is empty and `false` otherwise.
   */
  Status is_empty_bucket(const URI& uri, bool* is_empty) const;

  /**
   * Retrieves all the URIs that have the first input as parent.
   *
   * @param parent The target directory to list.
   * @param uris The URIs that are contained in the parent.
   * @return Status
   */
  Status ls(const URI& parent, std::vector<URI>* uris) const;

  /**
   * Retrieves all the entries contained in the parent.
   *
   * @param parent The target directory to list.
   * @return All entries that are contained in the parent
   */
  tuple<Status, optional<std::vector<filesystem::directory_entry>>>
  ls_with_sizes(const URI& parent) const;

  /**
   * Renames a file.
   *
   * @param old_uri The old URI.
   * @param new_uri The new URI.
   * @return Status
   */
  Status move_file(const URI& old_uri, const URI& new_uri);

  /**
   * Renames a directory.
   *
   * @param old_uri The old URI.
   * @param new_uri The new URI.
   * @return Status
   */
  Status move_dir(const URI& old_uri, const URI& new_uri);

  /**
   * Copies a file.
   *
   * @param old_uri The old URI.
   * @param new_uri The new URI.
   * @return Status
   */
  Status copy_file(const URI& old_uri, const URI& new_uri);

  /**
   * Copies directory.
   *
   * @param old_uri The old URI.
   * @param new_uri The new URI.
   * @return Status
   */
  Status copy_dir(const URI& old_uri, const URI& new_uri);

  /**
   * Reads from a file.
   *
   * @param uri The URI of the file.
   * @param offset The offset where the read begins.
   * @param buffer The buffer to read into.
   * @param nbytes Number of bytes to read.
   * @param use_read_ahead Whether to use the read-ahead cache.
   * @return Status
   */
  Status read(
      const URI& uri,
      uint64_t offset,
      void* buffer,
      uint64_t nbytes,
      bool use_read_ahead = true);

  /**
   * Reads multiple regions from a file.
   *
   * @param uri The URI of the file.
   * @param regions The list of regions to read. Each region is a tuple
   *    `(file_offset, dest_buffer, nbytes)`.
   * @param thread_pool Thread pool to execute async read tasks to.
   * @param tasks Vector to which new async read tasks are pushed.
   * @param use_read_ahead Whether to use the read-ahead cache.
   * @return Status
   */
  Status read_all(
      const URI& uri,
      const std::vector<tuple<uint64_t, Tile*, uint64_t>>& regions,
      ThreadPool* thread_pool,
      std::vector<ThreadPool::Task>* tasks,
      bool use_read_ahead = true);

  /**
   * Reads multiple regions from a file with no batching.
   *
   * @param uri The URI of the file.
   * @param regions The list of regions to read. Each region is a tuple
   *    `(file_offset, dest_buffer, nbytes)`.
   * @param thread_pool Thread pool to execute async read tasks to.
   * @param tasks Vector to which new async read tasks are pushed.
   * @param use_read_ahead Whether to use the read-ahead cache.
   * @return Status
   */
  Status read_all_no_batching(
      const URI& uri,
      const std::vector<tuple<uint64_t, Tile*, uint64_t>>& regions,
      ThreadPool* thread_pool,
      std::vector<ThreadPool::Task>* tasks,
      bool use_read_ahead = true);

  /** Checks if a given filesystem is supported. */
  bool supports_fs(Filesystem fs) const;

  /** Checks if the backend required to access the given URI is supported. */
  bool supports_uri_scheme(const URI& uri) const;

  /**
   * Syncs (flushes) a file. Note that for S3 this is a noop.
   *
   * @param uri The URI of the file.
   * @return Status
   */
  Status sync(const URI& uri);

  /**
   * Opens a file in a given mode.
   *
   *
   * @param uri The URI of the file.
   * @param mode The mode in which the file is opened:
   *     - READ <br>
   *       The file is opened for reading. An error is returned if the file
   *       does not exist.
   *     - WRITE <br>
   *       The file is opened for writing. If the file exists, it will be
   *       overwritten.
   *     - APPEND <b>
   *       The file is opened for writing. If the file exists, the write
   *       will start from the end of the file. Note that S3 does not
   *       support this operation and, thus, an error will be thrown in
   *       that case.
   * @return Status
   */
  Status open_file(const URI& uri, VFSMode mode);

  /**
   * Closes a file, flushing its contents to persistent storage.
   *
   * @param uri The URI of the file.
   * @return Status
   */
  Status close_file(const URI& uri);

  /**
   * Writes the contents of a buffer into a file.
   *
   * @param uri The URI of the file.
   * @param buffer The buffer to write from.
   * @param buffer_size The buffer size.
   * @return Status
   */
  Status write(const URI& uri, const void* buffer, uint64_t buffer_size);

  /**
   * Used in serialization to share the multipart upload state
   * among cloud executors during global order writes
   *
   * @param uri The file uri used as key in the internal map of the backend
   * @return A pair of status and VFS::MultiPartUploadState object.
   */
  std::pair<Status, std::optional<MultiPartUploadState>> multipart_upload_state(
      const URI& uri);

  /**
   * Used in serialization of global order writes to set the multipart upload
   * state in the internal maps of cloud backends during deserialization
   *
   * @param uri The file uri used as key in the internal map of the backend
   * @param state The multipart upload state info
   * @return Status
   */
  Status set_multipart_upload_state(
      const URI& uri, const MultiPartUploadState& state);

  /**
   * Used in remote global order writes to flush the internal
   * in-memory buffer for an URI that backends maintain to modulate the
   * frequency of multipart upload requests.
   *
   * @param uri The file uri identifying the backend file buffer
   * @return Status
   */
  Status flush_multipart_file_buffer(const URI& uri);

 private:
  /* ********************************* */
  /*        PRIVATE DATATYPES          */
  /* ********************************* */

  /**
   * Helper type holding information about a batched read operation.
   */
  struct BatchedRead {
    /** Construct a BatchedRead consisting of the single given region. */
    BatchedRead(const tuple<uint64_t, Tile*, uint64_t>& region) {
      offset = std::get<0>(region);
      nbytes = std::get<2>(region);
      regions.push_back(region);
    }

    /** Offset of the batch. */
    uint64_t offset;

    /** Number of bytes in the batch. */
    uint64_t nbytes;

    /**
     * Original regions making up the batch. Vector of tuples of the form
     * (offset, dest_buffer, nbytes).
     */
    std::vector<tuple<uint64_t, Tile*, uint64_t>> regions;
  };

  /**
   * Represents a sub-range of data within a URI file at a
   * specific file offset.
   */
  struct ReadAheadBuffer {
    /* ********************************* */
    /*            CONSTRUCTORS           */
    /* ********************************* */

    /** Value Constructor. */
    ReadAheadBuffer(const uint64_t offset, Buffer&& buffer)
        : offset_(offset)
        , buffer_(std::move(buffer)) {
    }

    /** Move Constructor. */
    ReadAheadBuffer(ReadAheadBuffer&& other)
        : offset_(other.offset_)
        , buffer_(std::move(other.buffer_)) {
    }

    /* ********************************* */
    /*             OPERATORS             */
    /* ********************************* */

    /** Move-Assign Operator. */
    ReadAheadBuffer& operator=(ReadAheadBuffer&& other) {
      offset_ = other.offset_;
      buffer_ = std::move(other.buffer_);
      return *this;
    }

    DISABLE_COPY_AND_COPY_ASSIGN(ReadAheadBuffer);

    /* ********************************* */
    /*             ATTRIBUTES            */
    /* ********************************* */

    /** The offset within the associated URI. */
    uint64_t offset_;

    /** The buffered data at `offset`. */
    Buffer buffer_;
  };

  /**
   * An LRU cache of `ReadAheadBuffer` objects keyed by a URI string.
   */
  class ReadAheadCache : public LRUCache<std::string, ReadAheadBuffer> {
   public:
    /* ********************************* */
    /*     CONSTRUCTORS & DESTRUCTORS    */
    /* ********************************* */

    /** Constructor. */
    ReadAheadCache(const uint64_t max_cached_buffers)
        : LRUCache(max_cached_buffers) {
    }

    /** Destructor. */
    virtual ~ReadAheadCache() = default;

    /* ********************************* */
    /*                API                */
    /* ********************************* */

    /**
     * Attempts to read a buffer from the cache.
     *
     * @param uri The URI associated with the buffer to cache.
     * @param offset The offset that buffer starts at within the URI.
     * @param buffer The buffer to cache.
     * @param nbytes The number of bytes within the buffer.
     * @param success True if `buffer` was read from the cache.
     * @return Status
     */
    Status read(
        const URI& uri,
        const uint64_t offset,
        void* const buffer,
        const uint64_t nbytes,
        bool* const success) {
      assert(success);
      *success = false;

      // Store the URI's string representation.
      const std::string uri_str = uri.to_string();

      // Protect access to the derived LRUCache routines.
      std::lock_guard<std::mutex> lg(lru_mtx_);

      // Check that a cached buffer exists for `uri`.
      if (!has_item(uri_str))
        return Status::Ok();

      // Store a reference to the cached buffer.
      const ReadAheadBuffer* const ra_buffer = get_item(uri_str);

      // Check that the read offset is not below the offset of
      // the cached buffer.
      if (offset < ra_buffer->offset_)
        return Status::Ok();

      // Calculate the offset within the cached buffer that corresponds
      // to the requested read offset.
      const uint64_t offset_in_buffer = offset - ra_buffer->offset_;

      // Check that both the start and end positions of the requested
      // read range reside within the cached buffer.
      if (offset_in_buffer + nbytes > ra_buffer->buffer_.size())
        return Status::Ok();

      // Copy the subrange of the cached buffer that satisfies the caller's
      // read request back into their output `buffer`.
      std::memcpy(
          buffer,
          static_cast<uint8_t*>(ra_buffer->buffer_.data()) + offset_in_buffer,
          nbytes);

      // Touch the item to make it the most recently used item.
      touch_item(uri_str);

      *success = true;
      return Status::Ok();
    }

    /**
     * Writes a cached buffer for the given uri.
     *
     * @param uri The URI associated with the buffer to cache.
     * @param offset The offset that buffer starts at within the URI.
     * @param buffer The buffer to cache.
     * @return Status
     */
    Status insert(const URI& uri, const uint64_t offset, Buffer&& buffer) {
      // Protect access to the derived LRUCache routines.
      std::lock_guard<std::mutex> lg(lru_mtx_);

      const uint64_t size = buffer.size();
      ReadAheadBuffer ra_buffer(offset, std::move(buffer));
      return LRUCache<std::string, ReadAheadBuffer>::insert(
          uri.to_string(), std::move(ra_buffer), size);
    }

   private:
    /* ********************************* */
    /*         PRIVATE ATTRIBUTES        */
    /* ********************************* */

    // Protects LRUCache routines.
    std::mutex lru_mtx_;
  };

  /* ********************************* */
  /*         PRIVATE ATTRIBUTES        */
  /* ********************************* */

#ifdef HAVE_AZURE
  Azure azure_;
#endif

#ifdef HAVE_GCS
  GCS gcs_;
#endif

#ifdef HAVE_S3
  S3 s3_;
#endif

#ifdef _WIN32
  Win win_;
#else
  Posix posix_;
#endif

#ifdef HAVE_HDFS
  tdb_unique_ptr<hdfs::HDFS> hdfs_;
#endif

  /** The class stats. */
  stats::Stats* stats_;

  /** The in-memory filesystem which is always supported */
  MemFilesystem memfs_;

  /**
   * Config.
   *
   * Note: This object is stored on the VFS for:
   * use of API 'tiledb_vfs_get_config'.
   * pass-by-reference initialization of filesystems' config_ member variables.
   **/
  Config config_;

  /** The set with the supported filesystems. */
  std::set<Filesystem> supported_fs_;

  /** Thread pool for compute-bound tasks. */
  ThreadPool* compute_tp_;

  /** Thread pool for io-bound tasks. */
  ThreadPool* io_tp_;

  /** Wrapper for tracking and canceling certain tasks on 'thread_pool' */
  CancelableTasks cancelable_tasks_;

  /** The read-ahead cache. */
  tdb_unique_ptr<ReadAheadCache> read_ahead_cache_;

  /* The VFS configuration parameters. */
  VFSParameters vfs_params_;

  /* ********************************* */
  /*          PRIVATE METHODS          */
  /* ********************************* */

  /**
   * Groups the given vector of regions to be read into a possibly smaller
   * vector of batched reads.
   *
   * @param regions Vector of individual regions to be read. Each region is a
   *    tuple `(file_offset, dest_buffer, nbytes)`.
   * @param batches Vector storing the batched read information.
   * @return Status
   */
  Status compute_read_batches(
      const std::vector<tuple<uint64_t, Tile*, uint64_t>>& regions,
      std::vector<BatchedRead>* batches) const;

  /**
   * Reads from a file by calling the specific backend read function.
   *
   * @param uri The URI of the file.
   * @param offset The offset where the read begins.
   * @param buffer The buffer to read into.
   * @param nbytes Number of bytes to read.
   * @param use_read_ahead Whether to use the read-ahead cache.
   * @return Status
   */
  Status read_impl(
      const URI& uri,
      uint64_t offset,
      void* buffer,
      uint64_t nbytes,
      bool use_read_ahead);

  /**
   * Executes a read, using the read-ahead cache as necessary.
   *
   * @param read_fn The read routine to execute.
   * @param uri The URI of the file.
   * @param offset The offset where the read begins.
   * @param buffer The buffer to read into.
   * @param nbytes Number of bytes to read.
   * @param use_read_ahead Whether to use the read-ahead cache.
   * @return Status
   */
  Status read_ahead_impl(
      const std::function<Status(
          const URI&, off_t, void*, uint64_t, uint64_t, uint64_t*)>& read_fn,
      const URI& uri,
      const uint64_t offset,
      void* const buffer,
      const uint64_t nbytes,
      const bool use_read_ahead);

  /**
   * Retrieves the backend-specific max number of parallel operations for VFS
   * read.
   */
  Status max_parallel_ops(const URI& uri, uint64_t* ops) const;
};

}  // namespace tiledb::sm

#endif  // TILEDB_VFS_H
