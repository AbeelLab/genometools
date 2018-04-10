package pelican

import atk.util.Tool

trait PelicanTool extends Tool{
  def getCCParams(cc: AnyRef) ={
    val m=(Map[String, Any]() /: cc.getClass.getDeclaredFields) { (a, f) =>
      f.setAccessible(true)
      a + (f.getName -> f.get(cc))
    }
    m.toList.sortBy(_._1)
  }
}