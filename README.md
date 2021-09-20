# processAndDemodulateEncodedDasData
processAndDemodulateEncodedDasData contains a collection of functions useful for extracting IRIG-B timecodes encoded in distributed acoustic sensor (DAS) data with a fiber stretcher.  
In particular, the top level function `processAndDemodulateEncodedDasData` will scan a differential phase directory for binary files, read and stack the dphase data associated with data marker channels, demodulate the AM signal, compute an IRIG-B bit string, and (finally) return the IRIG-B timestamp. 

This is done in several steps
1. Search Directory Tree for Binary Files (`dirRecursive3`)
2. Read blocks of lcm(SamplingRate, NumberOfPulses) / NumberOfPulses files. (e.g. for fs = 10000 Hz & nPulse = 16384 -> 625 binary files per block)  
  A. Binary data files are read with the C-MEX function `readPinnacleDasBinary`  
  B. A logical mask is applied to isolate the Data Marker Channels (Data Marker A & Data Marker B)  
  C. The data marker channels are stacked based on their average quality (e.g. IQ amplitude) as returned in the data quality block associated with the dPhase data (`stackSignal`)
3. An FIR filter is applied to the stacked data to remove low-frequency drifts (`fastFIRFilter`)
4. Product detectors are used to demodulate the AM signal (`MultiplyAndShiftDetector` and `HighBWLockInDetector` or their `gpuArray` equivalents `gpuMultiplyAndShiftDetector` and `gpuHighBWLockInDetector`)
5. Cross correlation is used to locate the seconds markers in the filtered signal (`xcorrFindSecondsMarkers`)
6. The demodulated AM signal is then split at the seconds markers and decoded to IRIG-B bit strings (`decodeIrigBSignalToBits` and `decodeIrigBSignalToBits4`)
7. The IRIG-B bit strings are padded and decoded to extract the timestamp (`allignIrigBBinaryWords` and `decodeIrigBSignalBitsFromArray`)
8. Timecode correction is done to backfill gaps (`TODO`)
9. File Time Stamps are Interpolated / Corrected based on the decoded IRIG-B timestamps
