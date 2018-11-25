using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data.SQLite;
using System.IO;

namespace ByonicDataParser
{
    class byonicDataParser
    {
        List<string> paths;
       
        string outputPath;

        public byonicDataParser(List<string> paths, string outputPath)
        {
            this.paths = paths;            
            this.outputPath = outputPath;
        }

        public void parse()
        {
            foreach(string path in paths)
            {
                onHighlightListItems(path);

                List<PSM> psms = new List<PSM>(); ;
                //Get Basic Information
                var browser = new SQLiteConnection(@"Data Source=" + path);
                browser.Open();

                var query = "SELECT ProteinsFoundPQMs.DebugText, PQMs.PeptideParse, PQMs.ProteinStartPosition, PQMs.PosteriorErrorProbability2, PQMs.PosteriorErrorProbability1,"+
                        "PQMs.FalseDiscoveryRate2, PQMs.FalseDiscoveryRate1, PQMs.FalseDiscoveryRateUnique2, PQMs.FalseDiscoveryRateUnique1, PQMs.PosteriorErrorProbability2_sum,"+
                        "PQMs.PosteriorErrorProbability1_sum, PQMs.CalcMz, ProteinsFoundPQMs.PQMsId, PQMsQueriesSummary.Intensity, Queries.ObservedMz, Queries.ScanTimeList, ProteinsFoundPQMs.ProteinRank, " +
                        "PQMs.Score, PQMs.DeltaScoreSeq, PQMs.DeltaScoreSeqMod, PQMs.Charge, PQMs.ObservedMh, PQMs.CalcMH, ProteinsFoundPQMs.Cleavage, Proteins.ProteinName,Proteins.Id, Queries.ScanNumberList "+
                        "FROM PQMs JOIN ProteinsFoundPQMs ON PQMS.Id = ProteinsFoundPQMs.PQMsId "+
                        "JOIN Queries ON PQMs.QueriesId = Queries.Id "+
                        "JOIN Proteins ON Proteins.Id = PQMs.ProteinsId " +
                        "JOIN PQMsQueriesSummary ON PQMsQueriesSummary.QueriesID = PQMs.QueriesID";
                var command = new SQLiteCommand(query, browser);
                var reader = command.ExecuteReader();

                while (reader.Read())
                {
                    //Values to get

                    //string PID = csv["PID"];                                                            //Ignored for now
                    string protRank = reader["ProteinRank"].ToString();
                    string pqmsID = reader["pqmsID"].ToString();
                    string sequence = reader["DebugText"].ToString();
                    string peptidesToBeParsed = reader["PeptideParse"].ToString();
                    int peptideStartPosition = int.Parse(reader["ProteinStartPosition"].ToString());   
                    double PEP2D = double.Parse(reader["PosteriorErrorProbability2"].ToString());
                    double PEP1D = double.Parse(reader["PosteriorErrorProbability1"].ToString());
                    double logProb = Math.Abs(Math.Log10(PEP1D));                                 //Ignored for now
                    double score = double.Parse(reader["Score"].ToString());
                    double deltaScore = double.Parse(reader["DeltaScoreSeq"].ToString());
                    double deltaModScore = double.Parse(reader["DeltaScoreSeqMod"].ToString());
                    int charge = int.Parse(reader["Charge"].ToString());
                    double mzObs = double.Parse(reader["ObservedMz"].ToString());
                    double mzCalc = double.Parse(reader["CalcMz"].ToString());                               
                    double obsMH = double.Parse(reader["ObservedMh"].ToString());
                    double calcMH = double.Parse(reader["CalcMH"].ToString());
                    string cleavage = reader["Cleavage"].ToString();                 
                    string proteinName = reader["ProteinName"].ToString();
                    int protID = int.Parse(reader["Id"].ToString());
                    string scanTime = reader["ScanTimeList"].ToString();                              

                    // For scan number
                    var scanNumberString = reader["ScanNumberList"].ToString();
                    var stringParts = scanNumberString.Split('=');
                    int scanNumber = int.Parse(stringParts[stringParts.Length-1]);
                    
                    double FDR2D = double.Parse(reader["FalseDiscoveryRate2"].ToString());
                    double FDR1D = double.Parse(reader["FalseDiscoveryRate1"].ToString());
                    double FDR2Dunique = double.Parse(reader["FalseDiscoveryRateUnique2"].ToString());
                    double FDR1Dunique = double.Parse(reader["FalseDiscoveryRateUnique1"].ToString());
                    double qvalue2D = double.Parse(reader["PosteriorErrorProbability2_sum"].ToString());
                    double qvalue1D = double.Parse(reader["PosteriorErrorProbability1_sum"].ToString());
                    double intensity = double.Parse(reader["intensity"].ToString());

                    PSM newPSM = new PSM(pqmsID, sequence, peptidesToBeParsed, peptideStartPosition, PEP2D, PEP1D, score, 
                        deltaScore, deltaModScore, charge, mzObs, mzCalc, obsMH, calcMH, cleavage, proteinName, protID, 
                        scanNumber, FDR2D, FDR1D, FDR2Dunique, FDR1Dunique, qvalue2D, qvalue1D, intensity, scanTime, protRank,
                        logProb);

                    psms.Add(newPSM);

                    //Console.ReadKey();
                }

                foreach(PSM psm in psms)
                {
                    query = "SELECT Modifications.Classification, Modifications.Composition, PQMsPeptideToModifications.ModificationsPeptidePosition,"+
                        "Modifications.AllowedSites, Modifications.MonoMassShiftTotal " +
                        "FROM PQMsPeptideToModifications "+ 
                        "JOIN PQMs ON PQMsPeptideToModifications.PQMsId = PQMs.Id "+
                        "JOIN Modifications ON Modifications.Id = PQMsPeptideToModifications.ModificationsId "+
                        "WHERE PQMsPeptideToModifications.PQMsId=" + psm.PQMsID;
                    command = new SQLiteCommand(query, browser);
                    reader = command.ExecuteReader();

                    string varMods = "";
                    string fixedMods = "";
                    string glycans = "";
                    string glycanPos = "";
                
                    while(reader.Read())
                    {
                        if (reader["Classification"].ToString().Equals("ngly") || reader["Classification"].ToString().Equals("ogly"))
                        {
                            glycans += reader["Composition"] + ";";
                            glycanPos += reader["ModificationsPeptidePosition"] + ";";

                            if (reader["Classification"].ToString().Equals("ngly"))
                            {
                                varMods += "N" + reader["ModificationsPeptidePosition"] + "(" + reader["AllowedSites"] + " / " 
                                            + reader["MonoMassShiftTotal"] + ");";
                            }
                            else
                            {
                                var position = Int32.Parse(reader["ModificationsPeptidePosition"].ToString());                                
                                varMods += psm.peptidesToBeParsed[position - 1] + reader["ModificationsPeptidePosition"].ToString() +
                                            "(" + reader["AllowedSites"] + " / " + reader["MonoMassShiftTotal"] + ");";
                            }
                                
                        }
                        else
                        {
                            if (reader["AllowedSites"].ToString().Equals("M"))
                            {
                                varMods += "M" + reader["ModificationsPeptidePosition"] + "(Oxidation / 15.9949);";
                            }

                            if(reader["AllowedSites"].ToString().Equals("NTerm E"))
                            {
                                varMods += "E" + reader["ModificationsPeptidePosition"] + "(Glu->pyro-Glu / -18.0106);";
                            }
                            if(reader["AllowedSites"].ToString().Equals("NTerm Q"))
                            {
                                varMods += "Q" + reader["ModificationsPeptidePosition"] + "(Gln->pyro-Glu / -17.026549);";
                            }
                            if (reader["AllowedSites"].ToString().Equals("N"))
                            {
                                varMods += "N" + reader["ModificationsPeptidePosition"] + "(Deamidated / 0.9840);";
                            }
                        }
                    }

                    
                    psm.modsFixed = fixedMods.TrimEnd(';');
                    psm.modsToBeParsed = varMods.TrimEnd(';');
                    psm.glycansToBeParsed = glycans.TrimEnd(';');
                    psm.glycanPositions = glycanPos.TrimEnd(';');

                }



                /**
                //TODO
                string modsToBeParsed = reader["Mods\r\n(variable)"].ToString();

                //TODO
                string glycansToBeParsed = reader["Glycans"].ToString();

                //TODO
                string glycanPositions = reader["Glycans\r\nPos."].ToString();

                //TODO
                string modsFixed = reader["Mods\r\n(fixed)"].ToString();
    **/

                string fileName = Path.GetFileNameWithoutExtension(path);

                StreamWriter writer = new StreamWriter(@outputPath + "\\" + fileName + "_exportedPeptides.txt");

                string headers = "ProtRank\tSequence\tPeptideParseFriendly\tPos.\tMods(variable)\tGlycans\tPEP2D\tPEP1D\t|Log Prob|\tScore\tDeltaScore\t" +
                                 "Delta Mod.Score\tz\tObs.m/z\tCalc.m/z\tppmerr.\tObs.MH\tCalc.MH\tCleavage\tGlycansPos.\tProteinName\tProt.Id\t" +
                                "ScanTime\tScan #\tMods(fixed)\tFDR2D\tFDR1D\tFDR uniq.2D\tFDR uniq.1D\tq-value2D\tq-value1D";

                writer.WriteLine(headers);

                foreach(PSM psm in psms)
                {
                    writer.WriteLine(psm.ToString());
                }

                writer.Close();
            
            }

            onFinish();
        }

        public EventHandler<HighlightEventArgs> highlightListItems;

        protected virtual void onHighlightListItems(string resultsFile)
        {
            var handler = highlightListItems;

            if(handler != null)
            {
                handler(this, new HighlightEventArgs(resultsFile));
            }
        }

        public EventHandler<EventArgs> finish;

        protected virtual void onFinish()
        {
            var handler = finish;

            if (handler != null)
            {
                handler(this, new EventArgs());
            }
        }
    }
}