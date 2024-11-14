#ifndef SELECTIONS_H
#define SELECTIONS_H

#include "ROOT/RDataFrame.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

RNode flagSelections(RNode df);
RNode triggerSelections(RNode df);
RNode electronSelections(RNode df);
RNode muonSelections(RNode df);
RNode leptonSelections(RNode df);
RNode higgsSelections(RNode df);
RNode WZSelections(RNode df);
RNode AK4Selections(RNode df);
RNode VBSJetsSelections(RNode df);

RNode applyPreSelections(RNode df);

#endif // SELECTIONS_H