using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScanAssigner
{
    class Glycan
    {
        private string _coreStructure;
        private List<SugarMoiety> _allSugars = new List<SugarMoiety>();
        private SugarConstants _sugarConstants = new SugarConstants();
        private Dictionary<string, List<SugarMoiety>> _allCombos;

        public Glycan(string CoreStructure)
        {
            _coreStructure = CoreStructure;
        }

        public void AddHexNAc(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.HexNAc);
            }
        }

        public void AddHex(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.Hex);
            }
        }

        public void AddFuc(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.Fuc);
            }
        }

        public void AddNeuAc(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.NeuAc);
            }
        }

        public void AddNeuGc(int Count)
        {
            for (int i = 0; i < Count; i++)
            {
                _allSugars.Add(_sugarConstants.NeuGc);
            }
        }

        public void GenerateCombinations()
        {
            var allCombos = GetAllCombos(_allSugars);
            GetUniqueCombos(allCombos);
            var t = "";
        }

        public string CoreStructure
        {
            get { return _coreStructure; }
        }

        public Dictionary<string, List<SugarMoiety>> AllFragments
        {
            get { return _allCombos; }
        }

        public void GetUniqueCombos(List<List<SugarMoiety>> allCombos)
        {
            _allCombos = new Dictionary<string, List<SugarMoiety>>();
            foreach (var list in allCombos)
            {
                Dictionary<string, int> names = new Dictionary<string, int>();
                foreach (var item in list)
                {
                    if (!names.ContainsKey(item.Name))
                    {
                        names.Add(item.Name, 0);
                    }
                    names[item.Name]++;
                }
                var nameKeys = names.Keys.ToList();
                nameKeys.OrderBy(x => x).ToList();
                var key = "";
                foreach (var name in nameKeys)
                {
                    key += string.Format(name + "({0})", names[name]);
                }
                if (!_allCombos.ContainsKey(key))
                {
                    _allCombos.Add(key, list);
                }
            }
        }

        public static List<List<T>> GetAllCombos<T>(List<T> list)
        {
            List<List<T>> result = new List<List<T>>();
            // head
            result.Add(new List<T>());
            result.Last().Add(list[0]);
            if (list.Count == 1)
                return result;
            // tail
            List<List<T>> tailCombos = GetAllCombos(list.Skip(1).ToList());
            tailCombos.ForEach(combo =>
            {
                result.Add(new List<T>(combo));
                combo.Add(list[0]);
                result.Add(new List<T>(combo));
            });
            return result;
        }
    }
}
