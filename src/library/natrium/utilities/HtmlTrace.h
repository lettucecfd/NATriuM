//============================================================================
// @file     HtmlTrace.h
// @author   Ottmar Kraemer-Fuhrmann <ottmar.kraemer-fuhrmann@scai.fraunhofer.de>
// @version  $Id$
// Copyright Copyright (C) 2013, Fraunhofer SCAI, Germany
// @brief    Stream integration test results to a html file.
//============================================================================

#ifndef HTML_TRACE_H
#define HTML_TRACE_H

#include <fstream>

/**
 * @short Stream integration test results to a html file.
 */
class HtmlTrace {

private:

	/// Stream for html output.
	std::ofstream m_html;

public:

	/// Constructor
	HtmlTrace();

	/// Destructor
	~HtmlTrace();

	/// Return stream.
	std::ofstream & getHtml() {
		return m_html;
	}
};

#endif // HTML_TRACE_H
