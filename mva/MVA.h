#ifndef MVA_H
#define MVA_H

#pragma once

#include <string>
#include <vector>

#include "argparser.hpp"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"

#include "TMVA/RInferenceUtils.hxx"
#include "TMVA/RReader.hxx"
#include "TMVA/RBDT.hxx"
#include "TMVA/RSofieReader.hxx"

using RNode = ROOT::RDF::RNode;

#endif // MVA_H