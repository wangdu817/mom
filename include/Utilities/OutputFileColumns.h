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

/**
 * @file OutputFileColumns.h
 * @brief Declaration of the OutputFileColumns utility class.
 */

#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace MOM
{

/**
 * @class OutputFileColumns
 * @brief Fixed-width, column-aligned ASCII output file writer for tabular data.
 *
 * Manages a single output file whose columns are declared once before any data
 * are written.  Each column is assigned a label and a floating-point precision;
 * the class computes a consistent field width that accommodates both the header
 * label and the widest scientific-notation value that the column can produce.
 *
 * The header row is written automatically by Complete() and every subsequent
 * row is produced by pairing a NewRow() call with successive `operator<<`
 * insertions — one per declared column, in order.
 *
 * @par Usage lifecycle
 * The object enforces a strict three-phase protocol:
 *
 * **Phase 1 — Declaration** (before any data):
 * @code
 *   MOM::OutputFileColumns out("results.txt");
 *   out.AddColumn("Time(s)",   6);   // label, significant digits
 *   out.AddColumn("Temp(K)",   4);
 *   out.AddColumn("fv",        8);
 *   out.Complete();                  // writes the header row; no more AddColumn after this
 * @endcode
 *
 * **Phase 2 — Data rows** (repeated for every time step / cell):
 * @code
 *   out.NewRow();              // terminate the previous row (first call is a no-op)
 *   out << time << T << fv;   // exactly n_columns insertions, in column order
 * @endcode
 *
 * **Phase 3 — Cleanup** (optional; destructor does not auto-close):
 * @code
 *   out.Close();               // flushes, closes the file, resets state
 * @endcode
 *
 * @par Column tag format
 * Each column is tagged as `"Label(N)"` where N is the 1-based column index,
 * e.g. `"Time(s)(1)"`, `"Temp(K)(2)"`.  This makes columns unambiguous when
 * labels are long or when the file is post-processed programmatically.
 *
 * @par Field-width calculation
 * The allocated field width for column i is
 * @code
 *   width[i] = max( tag_length + tab, precision + 9, minimum )
 * @endcode
 * where `tab = 4`, `minimum = 16`, and `precision + 9` covers the widest
 * possible scientific-notation token (`-X.XXX...e+YYY` including sign,
 * decimal point, exponent marker, exponent sign, and a three-digit exponent).
 *
 * @note This class is **not** thread-safe.  If multiple threads write
 *       diagnostics, each should own a separate instance.
 *
 * @note The file is opened in truncate mode (`std::ios::out`); any
 *       pre-existing file at the given path is silently overwritten.
 */
class OutputFileColumns
{
public:

    /**
     * @brief Default constructor.
     *
     * Initialises all internal state to defaults.  The object is in an
     * unopened state; Open() must be called before AddColumn().
     */
    OutputFileColumns();

    /**
     * @brief Constructor — opens the output file immediately.
     *
     * Equivalent to default-constructing and then calling Open(@p file_name).
     *
     * @param file_name Path to the output file.  The file is created or
     *                  truncated at construction time.
     */
    explicit OutputFileColumns(const std::filesystem::path& file_name);

    /**
     * @brief Constructor — opens the output file immediately (string overload).
     *
     * Equivalent to default-constructing and then calling Open(@p file_name).
     * Prefer the `std::filesystem::path` overload for new code; this overload
     * exists for compatibility with legacy string-based path APIs.
     *
     * @param file_name Path to the output file as a plain string.
     */
    explicit OutputFileColumns(const std::string& file_name);

    /**
     * @brief Open (or reopen) the output file.
     *
     * Opens the file at @p file_name for writing in truncation mode, sets the
     * stream to scientific notation with left-aligned fields.  Any previously
     * open file is closed first.
     *
     * @param file_name Path to the output file.
     *
     * @note Calling Open() on an already-open and partially-written file resets
     *       the internal state via Default() and discards any unflushed content.
     */
    void Open(const std::filesystem::path& file_name);

    /**
     * @brief Open (or reopen) the output file (string overload).
     *
     * Prefer the `std::filesystem::path` overload for new code.
     *
     * @param file_name Path to the output file as a plain string.
     */
    void Open(const std::string& file_name);

    /**
     * @brief Close the output file and reset all internal state.
     *
     * After this call the object is equivalent to a default-constructed instance
     * and may be reused by calling Open() followed by a fresh AddColumn() /
     * Complete() / data-row sequence.
     *
     * @note If the file is already closed, this function is a no-op with respect
     *       to the stream but still resets all bookkeeping state.
     */
    void Close();

    /**
     * @brief Declare one output column.
     *
     * Must be called after Open() and before Complete().  Each call appends one
     * column to the right of any previously declared columns.
     *
     * The column tag written to the header is `"<label>(<N>)"` where N is the
     * 1-based column index assigned at this call.
     *
     * @param label     Human-readable column label (e.g. `"Temperature(K)"`).
     *                  May be any non-empty string; special characters are not
     *                  escaped.
     * @param precision Number of significant digits for floating-point values
     *                  written to this column.  Also determines the minimum
     *                  field width (see class-level field-width formula).
     *
     * @note Calling AddColumn() after Complete() is not detected here; the
     *       erroneous extra column will cause the iterator check in NewRow() or
     *       Complete() to throw on the next phase transition.
     */
    void AddColumn(const std::string& label, unsigned int precision);

    /**
     * @brief Finalise column declarations and write the header row.
     *
     * Locks the column layout (no further AddColumn() calls are permitted),
     * resets the column iterator to 0, and writes one header line to the file
     * containing all column tags at their computed widths.
     *
     * @throws std::runtime_error If Complete() has already been called on this
     *                            instance, or if no columns have been declared
     *                            via AddColumn().
     *
     * @note Call exactly once per file, after all AddColumn() calls and before
     *       the first NewRow().
     */
    void Complete();

    /**
     * @brief Terminate the current data row and advance to the next.
     *
     * Writes a newline to the file and resets the column iterator to 0.  The
     * very first call (when iterator_ == 0 and no values have been inserted)
     * is a no-op with respect to column-count validation — this allows callers
     * to unconditionally call NewRow() at the top of a loop without special-
     * casing the first iteration.
     *
     * @throws std::runtime_error If the number of values inserted since the
     *                            last NewRow() does not equal the number of
     *                            declared columns (i.e. operator<< was called
     *                            too few or too many times).
     *
     * @note A trailing NewRow() after the last data row is recommended: it
     *       ensures the final row is properly terminated in the file.
     */
    void NewRow();

    /**
     * @brief Insert one value into the current column slot and advance.
     *
     * Formats @p number using the precision and field width of the current
     * column (determined by the insertion order) and writes it to the stream.
     * Returns `*this` to allow chaining:
     *
     * @code
     *   out << time << temperature << volume_fraction;
     * @endcode
     *
     * @tparam T Arithmetic type.  Supported instantiations are explicitly
     *           provided in `OutputFileColumns.cpp`: `int`, `unsigned int`,
     *           `long`, `unsigned long`, `long long`, `unsigned long long`,
     *           `float`, `double`, and `long double`.
     *
     * @param  number The value to write.
     * @return Reference to `*this` for operator chaining.
     *
     * @warning No bounds check is performed on the column iterator.  Inserting
     *          more values than declared columns invokes undefined behaviour
     *          (out-of-bounds vector access).  Always pair operator<< calls
     *          with a preceding NewRow() and a subsequent NewRow() to let the
     *          runtime check fire at the row boundary.
     */
    template <typename T>
    OutputFileColumns& operator<<(const T& number);

private:

    /**
     * @brief Reset all members to their constructed defaults.
     *
     * Called by constructors, Open(), and Close() to ensure the object always
     * starts each file in a clean, consistent state.  Clears all column
     * metadata vectors and resets all counters and flags.
     */
    void Default();

private:

    std::ofstream fOut_; //!< Output file stream (scientific, left-aligned).

    unsigned int n_;        //!< Total number of declared columns (set by Complete()).
    unsigned int iterator_; //!< Index of the next column to be written [0, n_).

    unsigned int tab_;      //!< Minimum gap (characters) between a column tag and the next field start.
    unsigned int minimum_;  //!< Absolute minimum field width applied to every column regardless of tag or precision.

    std::vector<std::string>  tag_;       //!< Column tags as written in the header, formatted as "Label(N)".
    std::vector<unsigned int> width_;     //!< Computed field widths per column [characters].
    std::vector<unsigned int> precision_; //!< Per-column floating-point significant-digit count.

    bool                  status_completed_; //!< True after Complete() has been called; guards against re-initialisation.
    std::filesystem::path file_name_;        //!< Path of the currently open file (used in error messages).
};

} // namespace MOM
