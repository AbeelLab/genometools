/**
 * %HEADER%
 */
package pelican;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYBarDataset;

import be.abeel.graphics.GraphicsFileExport;
import be.abeel.jfreechart.JFreeChartWrapper;
import be.abeel.util.FrequencyMap;

/**
 * Utility methods for frequency maps.
 * 
 * @author Thomas Abeel
 */
public class FMPlot {

	
	  public static void plot(FrequencyMap freq, String file) { plot(freq,
	  file, false); }
	  
	  /** Exports the frequency map as a plot in a png file.
	 */
	
	  public static void plot(FrequencyMap freq, String file, boolean
	  countNormalization) {
	  
	  plot(freq, file, countNormalization, 0, 0);
	  
	  }
	  
	  public static void plot(FrequencyMap freq, String file, boolean
	  countNormalization, int lower, int upper) { plot("Frequency map", freq,
	  file, countNormalization, lower, upper);
	  
	  }
	  
	  public static void plot(String title, FrequencyMap freq, String file,
	  boolean countNormalization, int lower, int upper) { List<FrequencyMap> l
	  = new ArrayList<FrequencyMap>(); l.add(freq); plot(title, l, file,
	  countNormalization, lower, upper, null, "bin", "count");
	  
	  }
	 
	// public static void plot(String title, FrequencyMap freq, String file,
	// boolean countNormalization, int lower, int upper, String[] labels,
	// String xAxis, String yAxis) {
	// List<FrequencyMap> l = new ArrayList<FrequencyMap>();
	// l.add(freq);
	// plot(title, l, file, countNormalization, lower,
	// upper,labels,xAxis,yAxis);
	// }

	public static void plot(String title, List<FrequencyMap> list, String file,
			boolean countNormalization, int lower, int upper, String[] labels,
			String xAxis, String yAxis) {

		boolean setted=(lower!=0&&upper!=0);
		
		System.out.println(setted);
		int min = setted? lower : Integer.MAX_VALUE;
		int max = setted? upper : Integer.MIN_VALUE;

		// XYSeries xy = new XYSeries("Frequency");
		// SimpleHistogramDataset xy=new SimpleHistogramDataset("Frequency");
		int count = 0;
		DefaultXYDataset set = new DefaultXYDataset();
		for (FrequencyMap freq : list) {

			double[][] data = new double[2][freq.keySet().size()];
			int index = 0;

			for (int i : freq.keySet()) {
				if (!setted) {
					if (i > max)
						max = i;
					if (i < min)
						min = i;
				}
				data[0][index] = i;
				data[1][index] = freq.get(i);
				// xy.add(i, this.get(i));
				index++;
				System.out.println(min+"\t"+max);
			}

			if (countNormalization) {
				int totalCount = freq.totalCount();
				for (int i = 0; i < data[1].length; i++)
					data[1][i] /= totalCount;
			}
			// XYSeriesCollection col = new XYSeriesCollection();

			set.addSeries(labels != null ? labels[count++] : "freq" + count++,
					data);
		}
		// col.addSeries(xy);
		JFreeChart chart = ChartFactory.createXYBarChart(title, xAxis, false,
				yAxis, new XYBarDataset(set, 1.0), PlotOrientation.VERTICAL,
				false, false, false);
		if (countNormalization) {
			chart.getXYPlot().getRangeAxis().setRange(0, 1);
		}
		System.out.println(min+"\t"+max);
	
		chart.getXYPlot().getDomainAxis().setRange(min, max);
		chart.getXYPlot().setBackgroundPaint(Color.WHITE);
		XYBarRenderer xy = (XYBarRenderer) chart.getXYPlot().getRenderer();
		xy.setShadowVisible(false);
		// xy.setDrawBarOutline(false);
		// xy.setDefaultBarPainter(new StandardXYBarPainter());
		xy.setBarPainter(new StandardXYBarPainter());
		xy.setGradientPaintTransformer(null);
		// xy.setDefaultBarPainter(new StandardXYBarPainter());
		// System.out.println(xy.getBaseFillPaint());
		// System.out.println(xy.getBaseShape());
		// System.out.println(xy.getBaseStroke());
		// val colors = List(Color.black,Color.red, new Color(0, 0, 0xCC),
		// Color.GREEN, new Color(0xFF, 0x99, 0x33) )
		// for (i <- 0 until chart.getXYPlot().getSeriesCount()) {
		xy.setSeriesPaint(0, Color.DARK_GRAY);
		//
		// }

		GraphicsFileExport.exportPDF(new JFreeChartWrapper(chart), file, 800,
				600);

	}

}
