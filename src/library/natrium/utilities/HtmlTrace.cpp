//============================================================================
// @file     HtmlTrace.cpp
// @author   Ottmar Kraemer-Fuhrmann <ottmar.kraemer-fuhrmann@scai.fraunhofer.de>
// @version  $Id$
// Copyright Copyright (C) 2013, Fraunhofer SCAI, Germany
// @brief    Stream integration test results to a html file.
//============================================================================

#include "HtmlTrace.h"

#include <time.h>
#include <unistd.h>
#include <mpi.h>

#include "deal.II/base/utilities.h"

//============================================================================
// Stream integration test results to a html file.
//============================================================================
HtmlTrace::HtmlTrace(bool active) :
		m_active(active) {

	if (m_active) {
		// Open html file
		m_html.open("natrium.html");
		m_html << "<html>\n\n";
		m_html << "<head>\n";
		m_html << "<title>NATriuM Integration Test</title>\n";
		m_html << "</head>\n\n";
		m_html << "<body>\n";

		// Header
		m_html << "<h2>NATriuM Integration Test</h2>\n";

		// Date & time
		time_t actTime = time(0);
		tm *now = localtime(&actTime);
		char timeBuffer[80];
		strftime(timeBuffer, 80, "%d.%m.%Y / Time: %H:%M:%S", now);
		m_html << "<h3>Date: " << timeBuffer;

		// Host
		char host[129];
		gethostname(host, 128);
		m_html << " / Host: " << host << "</h3>\n\n";

		// Table with header lines
		m_html << "<table border=\"1\">\n";
		/*	m_html << "  <tr>\n";
		 m_html << "    <th></th>\n";
		 m_html << "    <th colspan=\"3\" bgcolor=\"lightgray\">Testcase</th>\n";
		 m_html << "    <th rowspan=\"2\" bgcolor=\"lightgray\">Grid<br>Size</th>\n";
		 m_html << "    <th rowspan=\"2\" bgcolor=\"lightgray\">Phase<br>Selection</th>\n";
		 m_html << "    <th></th>\n";
		 m_html << "    <th rowspan=\"2\" bgcolor=\"lightgray\">Newton<br>Iterations</th>\n";
		 m_html << "    <th colspan=\"3\" bgcolor=\"lightgray\">Time(msec)</th>\n";
		 */
		m_html << "  </tr>\n";

		m_html << "  <tr bgcolor=\"lightgray\">\n";
		m_html << "    <th>#</th>\n";
		m_html << "    <th>Test</th>\n";
		m_html << "    <th>Quantity (K)</th>\n";
		m_html << "    <th>Expected</th>\n";
		m_html << "    <th>Outcome</th>\n";
		m_html << "    <th>Success</th>\n";
		m_html << "    <th>Time / sec</th>\n";
		m_html << "    <th>Details</th>\n";
		m_html << "    <th>Error details</th>\n";
		m_html << "  </tr>\n\n";
	} /* if (0 == dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))*/
} // HtmlTrace

//============================================================================
// Destructor.
//============================================================================
HtmlTrace::~HtmlTrace() {
	if (m_active) {

		// Finish table
		m_html << "</table>\n\n";

		// Footer
		m_html << "</body>\n";
		m_html << "</html>\n";

		// Close html file
		m_html.close();

	} /*if (0 == dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)) */
} // ~HtmlTrace

//============================================================================
// @short Add a row to the html table.
//============================================================================
void addHtmlRow(std::ofstream & html, size_t testCase) {
// Add row to html trace
		html << "  <tr align=\"center\">\n";
		html << "    <td bgcolor=\"lightgray\">" << testCase << "</td>\n";


} // addHtmlRow

//============================================================================
// @short Close a row of the html table.
//============================================================================
void closeHtmlRow(std::ofstream & html) {
// Close row of html trace
		html << "  </tr>\n";
		html.flush();
} // closeHtmlRow

