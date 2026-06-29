/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Lab <https://www.creckmodeling.polimi.it>              |
|   Department of Chemistry, Materials, and Chemical Engineering          |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of the OpenSMOKEpp library.                         |
|                                                                         |
|   Copyright (C) 2026 Alberto Cuoci.                                     |
|                                                                         |
|   OpenSMOKEpp is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKEpp is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKEpp. If not, see <https://www.gnu.org/licenses/>.  |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "Utilities/OutputFileColumns.h"

#include <stdexcept>

namespace MOM
{

// ============================================================================
// Constructors
// ============================================================================

OutputFileColumns::OutputFileColumns()
{
    Default();
}

OutputFileColumns::OutputFileColumns(const std::filesystem::path& file_name)
{
    Default();
    Open(file_name);
    file_name_ = file_name;
}

OutputFileColumns::OutputFileColumns(const std::string& file_name)
{
    Default();
    Open(file_name);
    file_name_ = file_name;
}

// ============================================================================
// Default — internal state reset
// ============================================================================

void OutputFileColumns::Default()
{
    n_                = 0;
    iterator_         = 0;
    minimum_          = 16; // minimum field width: wide enough for a typical scientific value
    tab_              = 4;  // gap between column tag end and the start of the data field
    status_completed_ = false;
    file_name_.clear();

    tag_.clear();
    width_.clear();
    precision_.clear();
}

// ============================================================================
// Open / Close
// ============================================================================

void OutputFileColumns::Open(const std::filesystem::path& file_name)
{
    file_name_ = file_name;

    // Truncate any existing file; set stream flags once for the lifetime of the file.
    fOut_.open(file_name, std::ios::out);
    fOut_.setf(std::ios::scientific); // values written as X.XXXe+YY
    fOut_.setf(std::ios::left);       // left-align within each field width
}

void OutputFileColumns::Open(const std::string& file_name)
{
    file_name_ = file_name;

    fOut_.open(file_name.c_str(), std::ios::out);
    fOut_.setf(std::ios::scientific);
    fOut_.setf(std::ios::left);
}

void OutputFileColumns::Close()
{
    if (fOut_.is_open())
        fOut_.close();

    // Reset bookkeeping so the object may be reused for a different file.
    Default();
}

// ============================================================================
// Column declaration
// ============================================================================

void OutputFileColumns::AddColumn(const std::string& label, const unsigned int precision)
{
    // Build the column tag as "Label(N)" where N is the 1-based column index.
    // The index is stringified via a stringstream to avoid sprintf/to_string
    // ambiguity on all supported platforms.
    iterator_++;
    std::stringstream number;
    number << iterator_;
    const std::string tag = label + "(" + number.str() + ")";
    tag_.push_back(tag);

    // Field-width calculation:
    //   label_width    = tag character count + tab  (keeps column tag readable)
    //   precision_width = precision + 9              (covers the widest scientific token:
    //                                                 sign[1] + digit[1] + '.'[1] +
    //                                                 fractional[precision] + 'e'[1] +
    //                                                 exp-sign[1] + exp-digits[3] = precision+8,
    //                                                 plus one space guard = precision+9)
    //   width[i]       = max(label_width, precision_width, minimum_)
    const unsigned int label_width     = static_cast<unsigned int>(tag.size()) + tab_;
    const unsigned int precision_width = precision + 9u;
    width_.push_back(std::max(label_width, std::max(precision_width, minimum_)));
    precision_.push_back(precision);
}

void OutputFileColumns::Complete()
{
    if (status_completed_)
    {
        std::stringstream message;
        message << "Filename: " << file_name_ << '\n'
                << "OutputFileColumns: Complete() has already been called; "
                   "columns cannot be re-declared.";
        throw std::runtime_error(message.str());
    }

    if (iterator_ == 0)
    {
        std::stringstream message;
        message << "Filename: " << file_name_ << '\n'
                << "OutputFileColumns: no columns have been declared via AddColumn().";
        throw std::runtime_error(message.str());
    }

    // Lock the layout and write the header row.
    status_completed_ = true;
    n_                = iterator_;
    iterator_         = 0;

    for (unsigned int i = 0; i < n_; i++)
        fOut_ << std::setprecision(precision_[i]) << std::setw(width_[i]) << tag_[i];

    fOut_ << std::endl;
}

// ============================================================================
// Data row interface
// ============================================================================

void OutputFileColumns::NewRow()
{
    // Allow the very first call (iterator_ == 0, no values written yet) to pass
    // through silently so callers can unconditionally place NewRow() at the top
    // of a data loop without special-casing the first iteration.
    if (iterator_ != n_ && iterator_ != 0)
    {
        std::stringstream message;
        message << "Filename: " << file_name_ << '\n'
                << "OutputFileColumns: row closed with " << iterator_
                << " value(s) but " << n_ << " column(s) were declared.";
        throw std::runtime_error(message.str());
    }
    fOut_ << std::endl;
    iterator_ = 0;
}

// ============================================================================
// operator<< — value insertion
// ============================================================================

template <typename T>
OutputFileColumns& OutputFileColumns::operator<<(const T& number)
{
    // Write the value formatted to the precision and width of the current column,
    // then advance the iterator so the next insertion targets the next column.
    fOut_ << std::setprecision(precision_[iterator_]) << std::setw(width_[iterator_]) << number;
    iterator_++;
    return *this; // return *this enables chaining: out << a << b << c;
}

// ============================================================================
// Explicit template instantiations
// ============================================================================
//
// operator<< is defined in this .cpp, not in the header, so the compiler
// cannot generate instantiations on demand at call sites.  The explicit
// instantiations below make all commonly used numeric types available to
// translation units that include only the header.
//
// To support an additional type T, add:
//   template OutputFileColumns& OutputFileColumns::operator<<(const T&);

// Integer types
template OutputFileColumns& OutputFileColumns::operator<<(const int&);
template OutputFileColumns& OutputFileColumns::operator<<(const unsigned int&);
template OutputFileColumns& OutputFileColumns::operator<<(const long&);
template OutputFileColumns& OutputFileColumns::operator<<(const unsigned long&);
template OutputFileColumns& OutputFileColumns::operator<<(const long long&);
template OutputFileColumns& OutputFileColumns::operator<<(const unsigned long long&);

// Floating-point types
template OutputFileColumns& OutputFileColumns::operator<<(const float&);
template OutputFileColumns& OutputFileColumns::operator<<(const double&);
template OutputFileColumns& OutputFileColumns::operator<<(const long double&);

} // namespace MOM
