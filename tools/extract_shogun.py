"""
Extract SVM parameters from the original Shogun 2.1.0 model pickles.

THIS SCRIPT RUNS UNDER PYTHON 2.7 WITH THE LEGACY SHOGUN BUILD. It is kept for
provenance and reproducibility only; nothing in the application needs it at
runtime. Its output (models/raw/*.npz) is committed to this repository, so the
legacy environment is not required to build or run SplashRNA v2.

Usage on the old server:

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/.../softwares/lib
    export PYTHONPATH=$HOME/.../softwares/splashrna/:$HOME/.../lib/python2.7/dist-packages
    python extract_shogun.py <model.bz2> <out.npz>

What it pulls out:
    alphas      support vector coefficients
    sv_idx      indices of the support vectors within the training set
    bias        decision function offset
    sv_strings  the 22-mer string of each support vector
    subkernel_weights

The kernel hyper-parameters (degree 10, shift 1, spectrum order 3, ...) are NOT
written here; they were read from a full ASCII serialization of the model and are
recorded as constants in tools/reference.py. See docs/MODEL-ARCHAEOLOGY.md.
"""

import sys

import numpy

# The pickles reference __main__.ShRNAtor, so the class must be importable there.
from sh_rna import ShRNAtor  # noqa: F401  (registered below)
import helper
import __main__

__main__.ShRNAtor = ShRNAtor


def extract(path, out_path):
    model = helper.load(path)
    svm = model.svm

    alphas = numpy.asarray(svm.get_alphas(), dtype=numpy.float64)
    sv_idx = numpy.asarray(svm.get_support_vectors(), dtype=numpy.int64)
    bias = float(svm.get_bias())

    import modshogun as sg

    combined = sg.CombinedKernel.obtain_from_generic(svm.get_kernel())
    features = combined.get_lhs()
    string_features = features.get_first_feature_obj()

    strings = []
    for i in range(string_features.get_num_vectors()):
        vec = numpy.asarray(string_features.get_feature_vector(i))
        if vec.dtype.kind in "iu":
            strings.append("".join(chr(int(c)) for c in vec))
        else:
            strings.append("".join(vec.astype("S1").tolist()))

    sv_strings = [strings[i] for i in sv_idx]

    numpy.savez_compressed(
        out_path,
        alphas=alphas,
        sv_idx=sv_idx,
        bias=numpy.array([bias]),
        sv_strings=numpy.array(sv_strings),
        subkernel_weights=numpy.asarray(
            combined.get_subkernel_weights(), dtype=numpy.float64
        ),
    )
    print("wrote %s: %d support vectors, bias %.12f" % (out_path, len(sv_idx), bias))


if __name__ == "__main__":
    extract(sys.argv[1], sys.argv[2])
