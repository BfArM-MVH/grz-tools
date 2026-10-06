# S3 upload integrity

`S3MultipartUploader` checks each upload with the Content-MD5 header and the ETag.
Every S3 backend that the GRZs use supports these two mechanisms.
The uploader does not use the Checksum API, because several of these backends ignore it or reject parts of it.
A file smaller than one part goes as one `PutObject` instead of a multipart upload, with the same checks.

## Mechanisms

| Mechanism                             | Computed by                                                                   | The server                                                                             | Backends        |
| ------------------------------------- | ----------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- | --------------- |
| ETag                                  | the server: the MD5 of the part bytes                                         | stores it with the object, but compares it with nothing. The client has to compare it. | all             |
| Content-MD5 header                    | the client: the MD5 of the part bytes                                         | rejects a part whose MD5 differs, and does not store the value.                        | all             |
| Checksum API (`x-amz-checksum-<alg>`) | the client: a CRC or SHA checksum of the part bytes                           | rejects a part whose checksum differs, and stores the value with the object.           | some, see below |
| COMPOSITE or FULL_OBJECT checksum     | the client: a checksum of the whole object, sent with CompleteMultipartUpload | checks the assembled object against it.                                                | some, see below |

The ETag of a multipart object is a hash of hashes: the MD5 of the concatenated part MD5s, followed by `-<number of parts>`.
A COMPOSITE checksum is built the same way from the part checksums.
A FULL_OBJECT checksum covers the bytes of the whole object, and exists only for CRC algorithms.

## Backend support for the Checksum API

- Ceph: support differs between deployments. Some ignore the checksum headers and accept wrong checksums. Others check each part, but ignore the checksums at CompleteMultipartUpload.
- Dell ECS: checks each part, but rejects every top-level checksum at CompleteMultipartUpload, whatever its value.
- IBM Cloud Object Storage: only in releases after 2025.
- SeaweedFS: unknown.

## What the uploader checks

1. Each part carries a Content-MD5 header, and so does the PUT of a file smaller than one part. The server rejects bytes that changed on the way.
2. The uploader compares the ETag that the server returns for each part or PUT with the local MD5.
3. After CompleteMultipartUpload, the uploader compares the object ETag with the hash of hashes of the local part MD5s.

The ETag comparison does not rely on the server rejecting anything.
A mismatch raises `UploadIntegrityError`.
A mismatch in step 3, or in the ETag of the PUT, also deletes the stored object.

## Memory

An upload holds about (threads + 1) x part size in memory for each file.
Each thread holds the part that it uploads, and one more part fills up from the stream.
The default `multipart_chunksize` is 256 MiB, so 4 threads hold about 1.25 GiB.
A file that would need more than 1000 parts gets larger parts, see `calculate_s3_part_size()`.

## Limits

The uploader expects MD5 ETags.
AWS S3 returns other ETags for objects that it encrypts with SSE-KMS or SSE-C.
So every upload to such a bucket fails with `UploadIntegrityError`.
