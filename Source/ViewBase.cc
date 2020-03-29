#include <iostream>

#include <TROOT.h>

#include "Garfield/Plotting.hh"
#include "Garfield/ViewBase.hh"

namespace Garfield {

ViewBase::ViewBase(const std::string& name) :
    m_className(name) { 

  plottingEngine.SetDefaultStyle();
}

ViewBase::~ViewBase() {
  if (!m_hasExternalCanvas && m_canvas) delete m_canvas;
}

void ViewBase::SetCanvas(TCanvas* c) {
  if (!c) return;
  if (!m_hasExternalCanvas && m_canvas) {
    delete m_canvas;
    m_canvas = nullptr;
  }
  m_canvas = c;
  m_hasExternalCanvas = true;
}

TCanvas* ViewBase::GetCanvas() {
  if (!m_canvas) {
    m_canvas = new TCanvas();
    if (m_hasExternalCanvas) m_hasExternalCanvas = false;
  }
  return m_canvas;
}

std::string ViewBase::FindUnusedFunctionName(const std::string& s) const {
  int idx = 0;
  std::string fname = s + "_0";
  while (gROOT->GetListOfFunctions()->FindObject(fname.c_str())) {
    ++idx;
    fname = s + "_" + std::to_string(idx);
  }
  return fname;
}

std::string ViewBase::FindUnusedHistogramName(const std::string& s) const {
  int idx = 0;
  std::string hname = s + "_0";
  while (gDirectory->GetList()->FindObject(hname.c_str())) {
    ++idx;
    hname = s + "_" + std::to_string(idx);
  }
  return hname;
}

}
