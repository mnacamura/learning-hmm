Revision history for Haskell package learning-hmm
===

## Version 0.3.1.2
- Default the limit of Baum-Welch iteration to 10000 (in `baumWelch'`)

## Version 0.3.1.1
- Bug fix release

## Version 0.3.1.0
- Add function `baumWelch'` that performs the Baum-Welch algorithm and returns
  a result locally maximizing its likelihood. This behaviour is different from
  that of `baumWelch`, which returns a list of intermediate results.

## Version 0.3.0.0
- Add `Learning.IOHMM` which represents a class of input-output HMM
- Delete `Learning.HMM.new`

## Version 0.2.1.0
- Performance improvements with employing 'hmatrix'
- `withEmission` now does iterative re-estimation based on the Viterbi path

## Version 0.2.0.0
- Remove dependency on the 'logfloat' package
- Performance improvements

## Version 0.1.1.0
- Add function `init` for random initialization
- Add function `simulate` for running a Markov process

## Version 0.1.0.0
- Original version
