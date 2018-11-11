using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using CSMSL;
using CSMSL.Spectral;
using CSMSL.Proteomics;
using CSMSL.IO.Thermo;
using LumenWorks.Framework.IO.Csv;

namespace LFQForm
{
    class LFQProcessor
    {
        public static List<LFPeptide> targetPeptides;

        public LFQProcessor(string rawFilePath, string peptidePath)
        {

            var rawFile = new ThermoRawFile(rawFilePath);
            rawFile.Open();

            List<int> msOneScanNumbers = new List<int>();
            foreach (var scane in rawFile)
            {
                if (scane.MsnOrder == 1)
                {
                    msOneScanNumbers.Add(scane.SpectrumNumber);
                }
            }

            ReadInTargets(peptidePath, rawFile);

            if (rawFile.GetMzAnalyzer(msOneScanNumbers[0]).Equals(MZAnalyzerType.Orbitrap))
            {
                var firstSpectrumRT = rawFile.GetRetentionTime(rawFile.FirstSpectrumNumber);
                var lastSpectrumRT = rawFile.GetRetentionTime(rawFile.LastSpectrumNumber);

                foreach (var pep in targetPeptides)
                {
                    // Initial peak boundaries
                    var startTime = Math.Max(pep.parentMS1Time - 2, firstSpectrumRT);
                    var stopTime = Math.Min(pep.parentMS1Time + 2, lastSpectrumRT);

                    pep.startLookupTime = startTime;
                    pep.stopLookupTime = stopTime;
                    pep.FirstScan = rawFile.GetSpectrumNumber(pep.startLookupTime);
                    pep.LastScan = rawFile.GetSpectrumNumber(pep.stopLookupTime);
                    //10 ppm maass error for peak picking
                    pep.lookupRange = DoubleRange.FromPPM(pep.UserMZ, 10);
                }

                foreach (var ms1 in msOneScanNumbers)
                {
                    var spectrum = rawFile.GetSpectrum(ms1);
                    GetXICs(spectrum, ms1, rawFile.GetRetentionTime(ms1));

                }
            }

        }

        public static void GetXICs(ThermoSpectrum currentSpectrum, int specNumber, double rt)
        {
            // Not sure what this is doing
            List<LFPeptide> donePeptides = targetPeptides.Where(x => x.LastScan < specNumber).ToList();
            foreach( var pep in donePeptides)
            {
                pep.doneBuildingXIC = true;
            }

            List<LFPeptide> currPeptides = targetPeptides.Where(x => x.FirstScan <= specNumber && x.LastScan >= specNumber).ToList();
            foreach (var pep in currPeptides)
            {
                List<ThermoMzPeak> outPeaks = new List<ThermoMzPeak>();
                if(currentSpectrum.TryGetPeaks(pep.lookupRange, out outPeaks))
                {
                    var matchPeak = GetClosestPeak(outPeaks, pep.UserMZ);
                    var newRTPeak = new RTPeak(matchPeak, rt);
                    pep.XICLibrary.Add(newRTPeak);
                }
                else
                {
                    var newRTPeak = new RTPeak(pep.UserMZ, 0, rt);
                    pep.XICLibrary.Add(newRTPeak);
                }
            }
        }

        public static ThermoMzPeak GetClosestPeak(List<ThermoMzPeak> peaks, double expMZ)
        {
            double diff = double.MaxValue;
            ThermoMzPeak returnPeak = null;
            foreach(var peak in peaks)
            {
                var currDiff = Math.Abs(peak.MZ - expMZ);
                if(currDiff < diff)
                {
                    diff = currDiff;
                    returnPeak = peak;
                }
            }
            return returnPeak;
        }

        public void ReadInTargets(string pepsPath, ThermoRawFile rawfile)
        {
            targetPeptides = new List<LFPeptide>();
            using(var csv = new CsvReader(new StreamReader(pepsPath), true))
            {
                int fieldCount = csv.FieldCount;
                string[] headers = csv.GetFieldHeaders();

                /**  Field Headers Needed in Scan Assigner Output!!! **/
                //"Spectrum number"
                //"Peptide"
                //"Charge"
                //"Precursor Theoretical m/z (Th)

                while (csv.ReadNextRecord())
                {
                    var pep = new LFPeptide
                    {
                        ms2ScanNumber = int.Parse(csv["Spectrum Number"]),
                        sequence = csv["Peptide"],
                        charge = int.Parse(csv["Charge"]),
                        TheoreticalMZ = double.Parse(csv["Precursor Theoretical m/z (Th)"])
                    };
                    // Get MS1 scan number
                    pep.parentMS1 = rawfile.GetParentSpectrumNumber(pep.ms2ScanNumber);
                    // Get retention time
                    pep.parentMS1Time = rawfile.GetRetentionTime(pep.parentMS1);
                    // Get precursor MZ
                    pep.UserMZ = rawfile.GetPrecursorMz(pep.ms2ScanNumber);

                    targetPeptides.Add(pep);
                }
            }
        }
    }
}
