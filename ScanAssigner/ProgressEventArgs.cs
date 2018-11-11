using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ScanAssigner
{
    class ProgressEventArgs : EventArgs
    {
        public double Progress { get; private set; }

        public ProgressEventArgs(double progress)
        {
            Progress = progress;
        }
    }
}
