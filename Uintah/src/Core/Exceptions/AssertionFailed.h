/*
 * The MIT License
 *
 * Copyright (c) 1997-2018 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */



/*
 *  AssertionFailed.h: Exception for a failed assertion.  Note - this
 *    version takes only a char*.  There is a FancyAssertionFailed that
 *    takes std::string.  This is done to prevent include file pollution
 *    with std::string
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   July 1999
 *
 */

#ifndef Core_Exceptions_AssertionFailed_h
#define Core_Exceptions_AssertionFailed_h

#include <Core/Exceptions/Exception.h>

#include   <string>

namespace Uintah {
class AssertionFailed : public Exception {
public:
  AssertionFailed(const char* msg,
		  const char* file,
		  int line);
  AssertionFailed(const AssertionFailed&);
  virtual ~AssertionFailed();
  virtual const char* message() const;
  virtual const char* type() const;
protected:
private:
  std::string message_;
  AssertionFailed& operator=(const AssertionFailed&);
};
} // End namespace Uintah

#endif


