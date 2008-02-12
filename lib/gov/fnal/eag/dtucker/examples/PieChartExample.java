package gov.fnal.eag.dtucker.examples;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.data.general.DefaultPieDataset;
import java.io.File;
public class PieChartExample {
    public static void main(String[] args) {
        // Create a simple pie chart
        DefaultPieDataset pieDataset = new DefaultPieDataset();
        pieDataset.setValue("A", new Integer(75));
        pieDataset.setValue("B", new Integer(10));
        pieDataset.setValue("C", new Integer(10));
        pieDataset.setValue("D", new Integer(5));
        //JFreeChart chart = ChartFactory.createPieChart(
       JFreeChart chart = ChartFactory.createPieChart3D(
               "CSC408 Mark Distribution",
                pieDataset, 
                true, 
                true, 
                false);
        try {
            ChartUtilities.saveChartAsJPEG(new File("piechart.jpg"), chart, 500,
                300);
        } catch (Exception e) {
            System.out.println("Problem occurred creating chart.");
        }
    }
}