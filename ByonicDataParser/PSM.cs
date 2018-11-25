using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ByonicDataParser
{
    class PSM
    {

        //tring PID;                                                                      
        public string protRank;                                                                
        public string PQMsID;
        public string sequence;
        public string peptidesToBeParsed;
        public int peptideStartPosition;

        //TODO
        public string modsToBeParsed;

        //TODO
        public string glycansToBeParsed;

        public double PEP2D;
        public double PEP1D;
        double logProb;                                                                 
        public double score;
        public double deltaScore;
        public double deltaModScore;
        public int charge;
        public double mzObs;
        public double mzCalc;
        double ppmError;                                                                
        public double obsMH;
        public double calcMH;
        public string cleavage;

        //TODO
        public string glycanPositions;

        public string proteinName;
        public int protID;
        int scanNumber;

        public double intensity;
        public string scanTime;

        //TODO
        public string modsFixed;


        public double FDR2D;
        public double FDR1D;
        public double FDR2Dunique;
        public double FDR1Dunique;
        public double qvalue2D;
        public double qvalue1D;


        public PSM(string PQMsID, string sequence, string peptidesToBeParsed, int peptideStartPosition, 
                    double PEP2D, double PEP1D, double score, double deltaScore, double deltaModScore,
                    int charge, double mzObs, double mzCalc, double obsMH, double calcMH, string cleavage,
                    string proteinName, int protID, int scanNumber, double FDR2D, double FDR1D, double FDR2Dunique,
                    double FDR1Dunique, double qvalue2D, double qvalue1D, double intensity, string scanTime, string protRank, double logProb)
        {
            this.PQMsID = PQMsID;
            this.sequence = sequence;
            this.peptidesToBeParsed = peptidesToBeParsed;
            this.peptideStartPosition = peptideStartPosition;
            this.PEP2D = PEP2D;
            this.PEP1D = PEP1D;
            this.score = score;
            this.deltaScore = deltaScore;
            this.deltaModScore = deltaModScore;
            this.charge = charge;
            this.mzObs = mzObs;
            this.mzCalc = mzCalc;
            this.obsMH = obsMH;
            this.calcMH = calcMH;
            this.cleavage = cleavage;
            this.proteinName = proteinName;
            this.protID = protID;
            this.scanNumber = scanNumber;
            this.FDR1D = FDR1D;
            this.FDR2D = FDR2D;
            this.FDR1D = FDR1D;
            this.FDR2Dunique = FDR2Dunique;
            this.FDR1Dunique = FDR1Dunique;
            this.qvalue2D = qvalue2D;
            this.qvalue1D = qvalue1D;
            this.intensity = intensity;
            this.scanTime = scanTime;
            this.ppmError = Math.Abs(mzCalc - mzObs)/mzCalc * 1000000;
            this.protRank = protRank;
            this.logProb = logProb;
        }

        public override string ToString()
        {
            string returnString =  protRank + "\t" + sequence + "\t" + peptidesToBeParsed + "\t" + 
                                   peptideStartPosition + "\t" + modsToBeParsed + "\t" + glycansToBeParsed + 
                                   "\t" + PEP2D + "\t" + PEP1D + "\t" + logProb + "\t" + score + "\t" + deltaScore + "\t" + 
                                   deltaModScore + "\t" + charge + "\t" + mzObs + "\t" + mzCalc + "\t" + ppmError + 
                                   "\t" + obsMH + "\t" + calcMH + "\t" + cleavage + "\t" + glycanPositions + "\t" + 
                                   proteinName + "\t" + protID + "\t" + Math.Round(double.Parse(scanTime)/60,2)  + "\t" + 
                                   scanNumber + "\t" + modsFixed + "\t" + FDR2D + "\t" + FDR1D + "\t" + FDR2Dunique + 
                                   "\t" + FDR1Dunique + "\t" + qvalue2D + "\t" + qvalue1D ;

            return returnString;
        }
    }
}

