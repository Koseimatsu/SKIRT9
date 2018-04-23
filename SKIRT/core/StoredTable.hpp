/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef STOREDTABLE_HPP
#define STOREDTABLE_HPP

#include "Array.hpp"
#include "StoredTableImpl.hpp"
#include <array>

////////////////////////////////////////////////////////////////////

/** An instance of the StoredTable<N> class template provides access to the contents of a
    particular resource file in the SKIRT stored table format (i.e. a "stored table").

    Stored table file format
    ------------------------

    A stored table includes the names of the quantities on the axes (e.g. wavelength and grain
    size) and those being tabulated (e.g. absorption and scattering efficiencies), in addition to
    the grid points for each axis and the tabulated data values. All grid points and values are
    stored as binary data in the form of 64-bit floating-point numbers, and are always given in SI
    units. The format is designed so that it is easy to calculate the offset, relative to the start
    of the file, to any particular data value. More specifically, a stored table file is
    essentially a sequence of 8-byte data items. A data item can have one of three types:
        - string: 1 to 8 printable and non-whitespace 7-bit ASCII characters, padded with spaces
          to fill 8 bytes if needed;
        - unsigned integer: 64-bit integer in little-endian byte order;
        - floating point: 64-bit double (IEEE 754) in little-endian byte order.

    The overall layout is as follows:
        - SKIRT name/version tag
        - Endianness tag
        - numAxes
        - axisName (x numAxes)
        - axisUnit (x numAxes)
        - axisScale (x numAxes)
        - [ numPoints  axisPoint (x numPoints) ] (x numAxes)
        - numQuantities
        - quantityName (x numQuantities)
        - quantityUnit (x numQuantities)
        - quantityScale (x numQuantities)
        - value (x numQuantities x numPoints1 x ... x numPointsN)
        - end-of-file tag

    The values are ordered so that the quantity values for a particular point are next to each
    other, the first axis index varies most rapidly, and the last axis index varies least rapidly.

    The StoredTable<N> class template
    ---------------------------------

    The StoredTable<N> template parameter \em N specifies the number of axes in the stored table,
    and thus the number of axis values needed to retrieve a tabulated quantity from the table. Each
    StoredTable<N> instance represents a single tabulated quantity. Accessing the full contents of
    a stored table resource file with multiple tabulated quantities requires a seperate
    StoredTable<N> instance for each of those quantities.

    The default constructor creates an invalid stored table instance. The alternate constructor and
    the open() function associate a particular stored table resource file with the stored table
    instance. The number of axes in this stored table resource file must match the template
    parameter \em N. Also, the axis names and the corresponding units in the file, and one of the
    tabulated quantity names and its corresponding unit in the file, must match the information
    passed to the alternate constructor or the open() function. The destructor automatically
    releases the file association and any related resources.

    The parenthesis operator returns the quantity represented by the StoredTable<N> for the \em N
    specified axis values, interpolated from the tabulated values. Other functions offer access
    for specific purposes, such as constructing a cumulative distribution function along one axis,
    given values for the other axes.

    Implementation and performance
    ------------------------------

    A StoredTable<N> instance acquires a read-only memory map on the associated stored table
    resource file as opposed to actually reading the file contents into memory through regular file
    I/O operations. This has some important, mostly positive, consequences.

    Acquiring a memory map establishes a mapping between "pages" of system-defined size in the
    logical address space of a process and the contents of the "backing file", in this case the
    stored table resource file. This operation is simple and thus very fast. From then on, the
    operating system automatically loads pages from the backing file into physical memory as they
    become needed because the program addresses an item in the logical memory range of the page.
    Conversely, the operating system automatically removes pages from physical memory if available
    memory becomes tight. In effect, the operating system automatically manages a high-performance
    caching mechanism on stored tables.

    Three important use cases benefit greatly from this mechanism. Firstly, a large resource file
    can be left associated with a StoredTable<N> instance for the duration of the program, even if
    it is used only sporadically. When memory is tight, infrequently used portions of the data will
    automatically be removed from memory and reloaded later if needed. Secondly, there is little
    overhead in constructing a StoredTable<N> instance (and possibly destroying it shortly
    thereafter) even when the program needs only a small portion of the file contents. And thirdly,
    because all StoredTable<N> instances associated with a given stored table resource file share
    the same memory map on that file, using a seperate instance for each quantity in the stored
    table incurs very little overhead.

    Moreover, most operating systems share memory maps between processes. For a parallel program
    using MPI, this means that all processes running on the same compute node share a single memory
    copy of the resources they employ. Also, most operating systems keep the memory map caches
    alive between consecutive invocations of a program (assuming memory is available), increasing
    performance when, for example, interactively testing the program.

    On the downside, a program requesting a huge chunk of data from a large stored table in a
    serial fashion would run faster using regular file I/O, because the separate page loads take
    more time than sequentially reading data in bulk. More importantly, performance usually
    degrades rapidly (to the point where the program no longer performs any useful work) when the
    system is forced to constantly remove and reload pages because there is not enough memory to
    hold the data needed for a particular phase in the program. And finally, the run-time
    performance of a program becomes somewhat unpredicable because the speed of accessing resources
    depends heavily on the previous state of the operating system caches.
*/
template<size_t N> class StoredTable
{
    static_assert(N >= 1, "StoredTable number of axes must be at least 1");

    // ================== Constructing ==================

public:
    /** The default constructor constructs an invalid stored table instance. The user must call the
        open() function to associate the stored table instance with a particular stored table
        resource file. Calling any of the other functions before calling open() results in
        undefined behavior (usually a crash). */
    StoredTable() { }

    /** This alternate constructor constructs a stored table instance and immediately associates a
        given stored table resource file with it by calling the open() function. Refer to the
        open() function for a description of the arguments and of its operation. */
    StoredTable(const SimulationItem* item, string filename, string axes, string quantity)
    {
        open(item, filename, axes, quantity);
    }

    /** This function associates a given stored table resource file with the stored table instance.
        If such an association already exists, this function throws a fatal error. Conversely,
        calling any of the other functions before an association exists results in undefined
        behavior (usually a crash).

        The \em item argument specifies a simulation item in the hierarchy of the caller (usually
        the caller itself) used to retrieve an appropriate logger.

        The \em filename argument specifies the filename of the resource, without any directory
        segments. The resource file must have the ".stab" filename extension, which will be added
        to the specified filename if needed.

        First of all, the number of axes in this stored table resource file must match the template
        parameter \em N. Furthermore, the axes names in the resource file and the corresponding
        units for each axis must match the information specified in the \em axes argument. Finally,
        one of the tabulated quantity names in the resource file and its corresponding unit must
        match the information specified in the \em quantity argument. For a stored table resource
        file with multiple tabulated quantities, the \em quantity argument at the same time
        determines which of these quantities will be associated with the stored table instance.

        The string passed to the \em axes argument must have the syntax
        "name1(unit1),...,nameN(unitN)". In other words, a unit string between parenthesis follows
        each axis name, and the specifications for different axes are separated by a comma. For
        example, "lambda(m),a(m)". Whitespace is not allowed. The string passed to the \em quantity
        argument must have a similar syntax, for a single name/unit combination. Examples include
        "Llambda(W/m)", "Qabs(1)", and "h(J/m3)".

        In summary, this function (1) locates the specified stored table resource file, (2)
        acquires a memory map on the file, (3) verifies that the stored table matches all
        requirements, and (4) stores relevant information in data members. If any of these steps
        fail, the function throws a fatal error. */
    void open(const SimulationItem* item, string filename, string axes, string quantity)
    {
        StoredTable_Impl::open(N, item, filename, axes, quantity,
                               _filePath, _axBeg.begin(), &_qtyBeg, _axLen.begin(), &_qtyStep,
                               _axLog.begin(), &_qtyLog);
    }

    /** The destructor breaks the association with a stored table resource file established by the
        alternate constructor or the open() function, if there is any. In practice, this simply
        means releasing the memory map on the associated file. */
    ~StoredTable() { StoredTable_Impl::close(_filePath); }

    // ================== Accessing values ==================

public:

    /** This function returns the value of the quantity represented by this stored table for the
        specified axes values, interpolated over the grid points of the actual tabulated values in
        all dimensions. The function uses linear or logarithmic interpolation for the axes and
        quantity values according to the flags specified in the stored table. Out-of-range axes
        values are automatically clamped to the corresponding outer grid point. */
    template <typename... Values, typename = std::enable_if_t<StoredTable_Impl::isNumericArgList<N, Values...>()>>
    double operator()(Values... values) const
    {
        std::array<double, N> axValues = {{ static_cast<double>(values)... }};
        return axValues[0];
    }

    /** For a one-dimensional table only, this function returns the value of the quantity
        represented by the stored table for the specified axes value, interpolated over the grid
        points of the actual tabulated values. The function uses linear or logarithmic
        interpolation for the axis and quantity values according to the flags specified in the
        stored table. Out-of-range axes values are automatically clamped to the corresponding outer
        grid point. */
    template <typename Value, typename = std::enable_if_t<N==1 && StoredTable_Impl::isNumericArgList<1, Value>()>>
    double operator[](Value value) const
    {
        // get the index of the upper border of the axis grid bin containing the specified axis value
        size_t right = std::lower_bound(_axBeg[0], _axBeg[0]+_axLen[0], value) - _axBeg[0];

        // if the value is beyond the grid borders, adjust both the bin border and the value
        if (right == 0)
        {
            right++;
            value = _axBeg[0][0];
        }
        else if (right == _axLen[0])
        {
            right--;
            value = _axBeg[0][right];
        }

        // get the index of the lower border of the axis grid bin containing the specified axis value
        size_t left = right-1;

        // get the values
        double x = value;
        double x1 = _axBeg[0][left];
        double x2 = _axBeg[0][right];
        double f1 = valueAtIndices(left);
        double f2 = valueAtIndices(right);

        // if requested, compute logarithm of coordinate values
        if (_axLog[0])
        {
            x  = log10(x);
            x1 = log10(x1);
            x2 = log10(x2);
        }

        // perform logarithmic interpolation of function value if requested and the bordering values are positive
        bool logf = _qtyLog && f1>0 && f2>0;

        // compute logarithm of function values if required
        if (logf)
        {
            f1 = log10(f1);
            f2 = log10(f2);
        }

        // perform the interpolation
        double fx = f1 + ((x-x1)/(x2-x1))*(f2-f1);

        // compute the inverse logarithm of the resulting function value if required
        if (logf) fx = pow(10,fx);

        return fx;
    }

    // ================== Accessing the raw data ==================

private:
    /** This function returns a copy of the value at the specified N indices. There is no range
        checking. Out-of-range index values cause unpredictable behavior. */
    template <typename... Indices, typename = std::enable_if_t<StoredTable_Impl::isIntegralArgList<N, Indices...>()>>
    double valueAtIndices(Indices... indices) const
    {
       return _qtyBeg[flattenedIndex(indices...)];
    }

    /** This template function returns the flattened index in the underlying data array for the
        specified N indices. */
    template <typename... Indices, typename = std::enable_if_t<StoredTable_Impl::isIntegralArgList<N, Indices...>()>>
    size_t flattenedIndex(Indices... indices) const
    {
        std::array<size_t, N> indexes = {{ static_cast<size_t>(indices)... }};
        size_t result = indexes[N-1];
        for (size_t k = N-2; k<N; --k) result = result * _axLen[k] + indexes[k];
        return result * _qtyStep;
    }

    // ================== Data members ==================

private:
    string _filePath;                       // the canonical path to the associated stored table file
    std::array<const double*, N> _axBeg;    // pointer to first grid point for each axis
    const double* _qtyBeg;                  // pointer to first quantity value
    std::array<size_t, N> _axLen;           // number of grid points for each axis
    size_t _qtyStep;                        // step size from one quantity value to the next (1=adjacent)
    std::array<bool, N> _axLog;             // interpolation type (true=log, false=linear) for each axis
    bool _qtyLog;                           // interpolation type (true=log, false=linear) for quantity
};

////////////////////////////////////////////////////////////////////

#endif
