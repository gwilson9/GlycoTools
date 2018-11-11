using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ByonicDataParser
{
    class HighlightEventArgs : EventArgs
    {
        public string resultsFile { get; private set; }

        public HighlightEventArgs(string resultsFile)
        {            
            this.resultsFile = resultsFile;
        }
    }
}
