## Description
- `BCH_encoder.m` encodes the messages in `msg.txt` using the generator polynomial with the specified design properties and stores them in `codeword.txt`
- `BCH_decoder.m` decodes the received codewords in `rx.txt` which contain **erasures** and **errors** and estimates the original message and stores them in `decoderOutput.txt`
- `mul_poly.m` is MATLAB function file which is used during decoding

## Note
- This code has been written in **MATLAB R2016a** version
- It is not tested for other versions of MATLAB
