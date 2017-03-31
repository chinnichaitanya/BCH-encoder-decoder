## Description
- `BCH_encoder.m` is the program file for encoder
- `BCH_decoder.m` is the program file for decoder
- `degree_poly.m` is a MATLAB function file which returns the degree of the polynomial given as input in ascending order of terms with respect to the power
- `mul_poly.m` is MATLAB function file which multipies polynomials over Galois field
- `gf2exp.m` is a MATLAB function file which converts a **gf** format array to **exponential** format over the Galois field 
- `mx.txt` contains the message vectors one per each line
- `codeword.txt` contains the codewords one per each line corresponding to each message vector in `msg.txt`
- `rx.txt` contains received codewords with errors and erasures (denoted by **2**) one per each line
- `decoderOut_codeword.txt` contains the decoded/estimated codewords one per each line corresponding to each received codeword in `rx.txt`
- `decoderOut_msg.txt` contains the decoded/estimated messages one per each line corresponding to each received codeword in `rx.txt`
- `logfile.log` contains the intermediate table values of the simplified Berlekamp-Massey algorithm for each decoded codeword

## Instructions to run
### Encoding
- Add the message bit vectors to the file `msg.txt` with one message vector in each line
- Execute `BCH_encoder.m` and it will output the encoded message vectors to the file `codeword.txt` in the same order as the messages in `msg.txt`

### Decoding
- Add the received codewords to the file `rx.txt` with erasures represented with '2'
- Execute `BCH_decoder.m` and it will output the decoded/estimated codewords to the file `decoderOut_coderowd.txt` and the decoded/estimated message bit vectors to the file `decoderOut_msg.txt`
- It also prints the intermediate table values of the Berlekamp-Massey algorithm, for all the cases to `logfile.log`
- Note that if the decoder fails, the corresponding error message will be printed instead of the estimated codeword and the estimated message vector

## Note
- This code has been written in **MATLAB R2016a** version
- It is not tested for other versions of MATLAB
