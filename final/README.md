## Description
- `BCH_encoder.m` encodes the messages in `msg.txt` using the generator polynomial with the specified design properties and stores them in `codeword.txt`
- `BCH_decoder.m` decodes the received codewords in `rx.txt` which contain **erasures** and **errors** and estimates the original message and stores them in `decoderOutput.txt`
- The `decoderOutput.txt` file is updated such that the estimated codeword and the message vector are one after another for every received vecotr and it is updated with error-message if decoder couldn't decode the received vector
- `mul_poly.m` is MATLAB function file which is used during decoding to multiply polynomials over Galois field
- `gf2exp.m` is a MATLAB function file which converts a **gf** format array to **exponential** format over the Galois field 

## Instructions to run
# Encoding
- Add the message bit vectors to the file `msg.txt` with one message vector in each line
- Execute `BCH_encoder.m` and it will output the encoded message vectors to the file `codeword.txt` in the same order as the messages in `msg.txt`

# Decoding
- Add the received codewords to the file `rx.txt` with erasures represented with '2'
- Execute `BCH_decoder.m` and it will output the decoded/estimated codewords to the file `decoderOut_coderowd.txt` and the decoded/estimated message bit vectors to the file `decoderOut_msg.txt`
- It also prints the intermediate table values of the Berlekamp-Massey algorithm, for all the cases to `logfile.log`
- Note that if the decoder fails, the corresponding error message will be printed instead of the estimated codeword and the estimated message vector

## Note
- This code has been written in **MATLAB R2016a** version
- It is not tested for other versions of MATLAB
