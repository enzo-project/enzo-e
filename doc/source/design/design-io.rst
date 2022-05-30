*******************
Scalable I/O Design
*******************

In Cello, "scalable I/O" is currently implemented using the
``"output"`` Method.  For scalability, we want to be able to write
large amounts of data to disk files in an efficient manner, and in a way
that subsequently can be efficiently read back in.

Data to be written include field and particle data, both of which are
associated with blocks in the array of octree mesh.  For large
simulations, there may be millions of such blocks, which may
collectively contain terabytes of data.

Writing this data to disk files efficiently requires flexibility in
the design. Perhaps the most influential decision is how many disk
files to write.  Writing to fewer files simplifies managing the data,
and reduces the cost of opening and closing files; whereas writing to
more files can increase parallelism and make more efficient use of the
underlying parallel file system.

Our current design uses a user-specified subset of root-level blocks
for performing the writes.

