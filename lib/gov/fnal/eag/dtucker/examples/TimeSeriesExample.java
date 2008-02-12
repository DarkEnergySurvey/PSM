package gov.fnal.eag.dtucker.examples;

import java.io.File;
import java.io.IOException;
import org.jfree.*;
import org.jfree.chart.*;
import org.jfree.data.time.*;

public class TimeSeriesExample {
	public static void main(String[] args) {
		// Create a time series chart
		TimeSeries pop = new TimeSeries("Population", Day.class);
		pop.add(new Day(10, 1, 2004), 100);
		pop.add(new Day(10, 2, 2004), 150);
		pop.add(new Day(10, 3, 2004), 250);
		pop.add(new Day(10, 4, 2004), 275);
		pop.add(new Day(10, 5, 2004), 325);
		pop.add(new Day(10, 6, 2004), 425);
		TimeSeriesCollection dataset = new TimeSeriesCollection();
		dataset.addSeries(pop);
		JFreeChart chart = ChartFactory.createTimeSeriesChart(
				"Population of CSC408 Town",
				"Date",
				"Population",
				dataset,
				true,
				true,
				false);
		try {
			ChartUtilities.saveChartAsJPEG(new File("timeserieschart.jpg"), chart, 500, 300);
		} catch (IOException e) {
			System.err.println("Problem occurred creating chart.");
		}
	}
}
