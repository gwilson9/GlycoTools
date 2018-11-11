using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CSMSL.Chemistry;


namespace ScanAssigner
{
    class SugarMoiety
    {
        private string _name;
        private string _symbol;
        private ChemicalFormula _chemicalFormula;

        public SugarMoiety(string Name, string Symbol, string ChemicalFormula)
        {
            _name = Name;
            _symbol = Symbol;
            _chemicalFormula = new ChemicalFormula(ChemicalFormula);
        }

        public string Name
        {
            get { return _name; }
        }

        public string Symbol
        {
            get { return _symbol; }
        }

        public ChemicalFormula ChemicalFormula
        {
            get { return _chemicalFormula; }
        }
    }
}
