1. Basic tests : what they are used for
A basic test consists in a sequence one or more nucleotides have been added to. At those sites Augustus is asked to apply fs. 
The test succeeds if the applying Augustus with fs option enabled to the modified sequence produces the same gff you can obtain
by applying Augustus with fs option disabled at the original sequence.
This is intended to be a basic test to verify that fs at any site is performed without any side-effect to the preexisting code.

2. Basic test : how to create them
In order to insert artificial fs within any sequence, prepareTests can be used

3. Basic tests : how to run them
In order to run a test, overwrite the contents of bin/fs.hints with the site where fs is wanted to occurr and type in the following command line:
./augustus --fs=file --species=fly sampleFile.fasta > sampleFile.gff
The test succeeds if all features in the final gff have the same start/end as the original gff (except for +1 shift)
