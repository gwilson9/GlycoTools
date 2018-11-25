using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL.Spectral;
using CSMSL.IO.Thermo;
using CSMSL.Proteomics;

namespace ScanAssigner
{
    class PSM
    {
        public Peptide peptide { get; set; }
        public string PID { get; set; }
        public string protRank { get; set; }
        public string peptidesToBeParsed { get; set; }
        public int peptideStartPosition { get; set; }
        public string modsToBeParsed { get; set; }
        public string glycansToBeParsed { get; set; }
        public double PEP2D { get; set; }
        public double PEP1D { get; set; }
        public double logProb { get; set; }
        public double score { get; set; }
        public double deltaScore { get; set; }
        public double deltaModScore { get; set; }
        public int charge { get; set; }
        public double mzObs { get; set; }
        public double mzCalc { get; set; }
        public string ppmError { get; set; }
        public double obsMH { get; set; }
        public double calcMH { get; set; }
        public string cleavage { get; set; }
        public string glycanPositions { get; set; }
        public string proteinName { get; set; }
        public int protID { get; set; }
        public double scanTime { get; set; }
        public int scanNumber { get; set; }
        public string modsFixed { get; set; }
        public double FDR2D { get; set; }
        public double FDR1D { get; set; }
        public double FDR2Dunique { get; set; }
        public double FDR1Dunique { get; set; }
        public double qvalue2D { get; set; }
        public double qvalue1D { get; set; }
        public bool isGlycoPeptide { get; set; }
        public string fragmentation { get; set; }
        public double intensity;
        public string byonicSequence;
        public bool modsPassedCheck;
        public int masterScan;
        public bool positionPassedCheck;
        public string byonicIntensity;
        //public ThermoSpectrum spectrum { get; set; }

        public ThermoMzPeak Peak126;
        public ThermoMzPeak Peak138;
        public ThermoMzPeak Peak144;
        public ThermoMzPeak Peak168;
        public ThermoMzPeak Peak186;
        public ThermoMzPeak Peak204;
        public ThermoMzPeak Peak274;
        public ThermoMzPeak Peak292;
        public ThermoMzPeak Peak366;


        public PSM(Peptide peptide, string PID, string protRank, string peptidesToBeParsed, int peptideStartPosition, string modsToBeParsed, string glycansToBeParsed, double PEP2D, double PEP1D, double logProb, double score,
            double deltaScore, double deltaModScore, int charge, double mzObs, double mzCalc, string ppmError, double obsMH, double calcMH, string cleavage, string glycanPositions, string proteinName, int protID,
            double scanTime, int scanNumber, string modsFixed, double FDR2D, double FDR1D, double FDR2Dunique, double FDR1Dunique, double qvalue2D, double qvalue1D, bool isGlycoPeptide, string fragmentation, string byonicSequence, 
            bool modsPassedCheck, int masterScan, bool positionPassedCheck, string byonicIntensity)//, ThermoSpectrum spectrum)
        {
            this.peptide = peptide;
            this.PID = PID;
            this.protRank = protRank;
            this.peptidesToBeParsed = peptidesToBeParsed;
            this.peptideStartPosition = peptideStartPosition;
            this.modsToBeParsed = modsToBeParsed;
            this.glycansToBeParsed = glycansToBeParsed;
            this.PEP2D = PEP2D;
            this.PEP1D = PEP1D;
            this.logProb = logProb;
            this.score = score;
            this.deltaScore = deltaScore;
            this.deltaModScore = deltaModScore;
            this.charge = charge;
            this.mzObs = mzObs;
            this.mzCalc = mzCalc;
            this.ppmError = ppmError;
            this.obsMH = obsMH;
            this.calcMH = calcMH;
            this.cleavage = cleavage;
            this.glycanPositions = glycanPositions;
            this.proteinName = proteinName;
            this.protID = protID;
            this.scanTime = scanTime;
            this.scanNumber = scanNumber;
            this.modsFixed = modsFixed;
            this.FDR2D = FDR2D;
            this.FDR1D = FDR1D;
            this.FDR2Dunique = FDR2Dunique;
            this.FDR1Dunique = FDR1Dunique;
            this.qvalue2D = qvalue2D;
            this.qvalue1D = qvalue1D;
            this.isGlycoPeptide = isGlycoPeptide;
            this.fragmentation = fragmentation;
            this.byonicSequence = byonicSequence;
            this.modsPassedCheck = modsPassedCheck;
            this.masterScan = masterScan;
            this.byonicIntensity = byonicIntensity;
            this.positionPassedCheck = positionPassedCheck;
            //this.spectrum = spectrum;
        }
    }
}
