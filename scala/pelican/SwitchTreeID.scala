package pelican

import atk.util.Tool
import java.io.PrintWriter
import abeel.genometools.Main

object SwitchTreeID extends Main{

  override val version="""
    2014/11/28    Initial version
                  
    
    
    """
  override def main(args: Array[String]): Unit = {
    
    val inf="n:/Projects/TB-ARC/INDIA/jun9/jun9.rooted.nwk"
    
    val input=tLines(inf)
    
    assume(input.size==1)
    
    val map=tMap(tLines("n:/Projects/TB-ARC/INDIA/conversion.txt"),0,1,limitSplit=false)
    
    var mod=input(0)
    map.map(f=>{
      mod=mod.replaceAll(f._1, f._2)
      
      
    })
    val pw=new PrintWriter(inf+".renamed")
    pw.println(mod)
    pw.close
    
  }

}