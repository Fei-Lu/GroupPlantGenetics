package xiaohan.utils;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;

/**
 * java database utils
 *
 * @ author: yxh
 * @ created: 2021-08-06 : 9:41 AM
 */
public class DBUtils {
    private static final String URL = "jdbc:mysql://localhost:3306/demo_jdbc";
    private static final String NAME = "root";
    private static final String PASSWORD = "root";

    public static void main(String[] args) throws Exception {
        //1.加载驱动程序
        Class.forName("com.mysql.jdbc.Driver");
        //2.获得数据库的连接
        Connection conn = DriverManager.getConnection(URL, NAME, PASSWORD);
        //3.通过数据库的连接操作数据库，实现增删改查
        Statement stmt = conn.createStatement();
        ResultSet rs = stmt.executeQuery("select user_name,age from imooc_goddess");//选择import java.sql.ResultSet;
        while (rs.next()) {//如果对象中有数据，就会循环打印出来
            System.out.println(rs.getString("user_name") + "," + rs.getInt("age"));
        }
    }
}
